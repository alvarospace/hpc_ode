// CVODE INCLUDES
#include <app_magma/cvode_user.cuh>

// UTILS
#include <app_magma/utils.hpp>

// Include for memcpy
#include <cstring>


void cvode_run(const std::string& inputFile, const std::string& outputFile,
               const size_t pack_size) {

    /********** INPUT CONSTANTS ************/

    // Pressure Constant (Pa)
    realtype P = 101325.15;

    // Time step
    realtype dt = 1e-3;

    // Tolorences 
    realtype reltol = RCONST(1.0e-6);
    realtype abstol = RCONST(1.0e-10);
    

    /****************************/

    // TIMERS
    Utils::Timer simTime;

    /*********** READ CSV FILE and SET PACKAGE SIZE ************/
    std::shared_ptr<Utils::ThermData> mesh = Utils::readCsv(inputFile);
    mesh->info();
    
    // Number of points of the mesh
    size_t n_size = mesh->points;

    // Outer loop iterations
    size_t mod = n_size % BLOCKSIZE;
    size_t ext_it = n_size / BLOCKSIZE;
    // If the division is not exact, CPU will process the remainder
    

    std::cout << "BLOCKSIZE GPU: " << BLOCKSIZE << std::endl;
    std::cout << "Number of packages: " << ext_it << std::endl;
    std::cout << "Remainder points (will be process in serial in CPU): " << mod << std::endl;
    /************************************************************/

    // Check that number of species read is the same that in pyjac library
    if (mesh->nsp != NSP) {
        std::cout << "Error: Number of species in the input file does not match Pyjac species" << std::endl;
        exit(EXIT_FAILURE);
    }

    /********** CVODE INTEGRATION ************/

    /**** GPU ****/
    // Variables declarations
    int m_maxsteps = 10000;

    // // Problem size has to be "sunindextype", NSP * nBatch
    sunindextype nsp_GPU;

    /* Sundials context creation and automatic free memory (C++ feature) */
    sundials::Context sunctx;


    /* Cvode mem */
    void* cvode_mem = NULL;

    /* CUDA Memory helper (sundials help managing copies between CPU and device) */
    SUNMemoryHelper cuda_mem_help;
    cuda_mem_help = SUNMemoryHelper_Cuda(sunctx);

    /* Data variables */
    N_Vector y;
    realtype *yptr;
    SUNMatrix J;
    SUNLinearSolver LS;

    int retval;
    /*************/


    /**** CPU ****/
    N_Vector y_cpu;
    realtype *yptr_cpu;
    SUNMatrix J_cpu;
    SUNLinearSolver LS_cpu;
    void* cvode_mem_cpu = NULL;
    sundials::Context sunctx_cpu;
    int retval_cpu;

    /*************/
    simTime.tic();
    #pragma omp parallel num_threads(2){
        int id = omp_get_thread_num();

        // Integration GPU -> id thread = 0
        if (id == 0){
            for (int i = 0; i < ext_it; i++) {
                size_t gpu_points;

                // Limit the pack_size to the last iteration
                // because the number of points to process is the remainder
                if (i == (ext_it - 1) && mod > 0)
                    gpu_points = mod;
                else
                    gpu_points = pack_size;


                /* Start initialization rutine for GPU */

                // GPU PyJac memory allocation requirements
                mechanism_memory *h_mem, *d_mem;
                h_mem = new mechanism_memory;
                initialize_gpu_memory(gpu_points, &h_mem, &d_mem);
                

                // CVode User defined data to access in user supplied functions
                UserData *userData = new UserData;
                userData->Pressure = P;
                userData->h_mem = h_mem;
                userData->d_mem = d_mem;


                // Problem size for GPU
                nsp_GPU = NSP * gpu_points;
                userData->nEquations = nsp_GPU;
                userData->nSystems = gpu_points;

                // Vector allocation and initial conditions
                y = N_VNew_Cuda(nsp_GPU, sunctx);
                if (check_retval((void *)y, "N_VNew_CUDA", 0)) exit(EXIT_FAILURE);
                yptr = N_VGetArrayPointer(y);

                for (int j = 0; j < gpu_points; j++) {
                    // Index of the global point
                    size_t gIndex = i * pack_size + j;
                    size_t lIndex = j * NSP;

                    yptr[lIndex] = mesh->temp[gIndex];
                    std::memcpy(&yptr[lIndex] + 1, mesh->matSp[gIndex].data(), (NSP - 1) * sizeof(realtype));
                }
                
                // Copy initial data to GPU
                N_VCopyToDevice_Cuda(y);

                /* Call CVodeCreate to create the solver memory and specify the
                * Backward Differentiation Formula */
                cvode_mem = CVodeCreate(CV_BDF, sunctx);
                if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) exit(EXIT_FAILURE);


                /* CVode init dydt = dydt_cvode(t0, y) */
                retval = CVodeInit(cvode_mem, dydt_cvode, 0.0, y);
                if (check_retval(&retval, "CVodeInit", 1)) exit(EXIT_FAILURE);


                /* Save UserData structure internally in CVode */
                retval = CVodeSetUserData(cvode_mem, userData);
                if (check_retval(&retval, "CVodeSetUserData", 1)) exit(EXIT_FAILURE);

                /* Call CVodeSStolerances */
                retval = CVodeSStolerances(cvode_mem, reltol, abstol);
                if (check_retval(&retval, "CVodeSStolerances", 1)) exit(EXIT_FAILURE);


                /* Create SUNMatrix for use in linear solves */
                J = SUNMatrix_MagmaDenseBlock(gpu_points, NSP, NSP, SUNMEMTYPE_DEVICE, cuda_mem_help, NULL, sunctx);
                if(check_retval((void *)J, "SUNMatrix_MagmaDenseBlock", 0)) exit(EXIT_FAILURE);
                

                /* Create dense SUNLinearSolver (for magma library) */
                LS = SUNLinSol_MagmaDense(y, J, sunctx);
                if (check_retval((void *)LS, "SUNLinSol_MagmaDense", 0)) exit(EXIT_FAILURE);

                /* Attach the matrix and linear solver */
                retval = CVodeSetLinearSolver(cvode_mem, LS, J);
                if (check_retval(&retval, "CVodeSetLinearSolver", 1)) exit(EXIT_FAILURE);

                /* Set the user-supplied Jacobian routine Jac */
                retval = CVodeSetJacFn(cvode_mem, eval_jacob_cvode);
                if (check_retval(&retval, "CVodeSetJacFn", 1)) exit(EXIT_FAILURE);

                /* Max num steps of the method */
                retval = CVodeSetMaxNumSteps(cvode_mem, m_maxsteps);
                if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) exit(EXIT_FAILURE);


                // Run simulation
                realtype t0 = 0.0f;

                retval = CVode(cvode_mem, dt, y, &t0, CV_NORMAL);

                N_VCopyFromDevice_Cuda(y);

                /*
                *   - Calculation of the last element (mass convervation) of every point
                *   - Save data in mesh structure
                */
                std::vector<double> yLast(gpu_points, 0.0f);
                for (int j = 0; j < gpu_points; j++) {

                    // Index of the global point
                    size_t gIndex = i * pack_size + j;
                    size_t lIndex = j * NSP;

                    for (int k = 1; k < NSP; k++) {
                        yLast[lIndex] += yptr[lIndex + k];
                    }

                    mesh->temp[gIndex] = yptr[lIndex];
                    std::memcpy(mesh->matSp[gIndex].data(), &yptr[lIndex] + 1, (NSP-1) * sizeof(realtype));
                    mesh->matSp[gIndex][NSP-1] = yLast[lIndex];
                }


                /* Free Memory */
                free_gpu_memory(&h_mem, &d_mem);
                delete(h_mem);
                delete(userData);

                N_VDestroy(y);
                CVodeFree(&cvode_mem);
                SUNLinSolFree(LS);
                SUNMatDestroy(J);
            }
        } else if(mod > 0) {
            int first_cpu_it = ext_it * BLOCKSIZE;
            for (int i = first_cpu_it; i < first_cpu_it + mod; i++) {
                // TO DO Serial approach
            }
        }
    }
    simTime.toc();

    std::cout << "Calculation time: " << simTime.time() << " s" << std::endl;

    /* Write results for validation */
    Utils::writeCsv(mesh, outputFile); 
    std::cout << "AQUI LLEGO" << std::endl;
}

int main(int argc, char *argv[]) {
    std::string inputFile {"../../../ref_data/res_gri3.0.csv"}; 
    std::string outputFile {"results.csv"};
    size_t packSize {10};

     /* Command-line arguments logic */
    std::vector<std::string> args;
    if (argc > 1) {
        for (int i = 0; i < argc; i++) {
            args.push_back(std::string(argv[i]));
        }

        if (args[1].compare("--help") == 0 || args[1].compare("-h") == 0) {
            std::cout << "Usage: " << args[0] <<" [input file=../../../ref_data/res_gri3.0.csv] [output file=results.csv] [packSize=10]" << std::endl;
            return 0;
        } else {
            if (argc == 2) {
                inputFile = args[1];
            } else if (argc == 3) {
                inputFile = args[1];
                outputFile = args[2];
            } else if (argc == 4) {
                inputFile = args[1];
                outputFile = args[2];
                packSize = std::stod(args[3]);
            } else {
                std::cout << "Too many arguments..." << std::endl;
                return EXIT_FAILURE;
            }
        }
    }

    /* Run integrator */
    cvode_run(inputFile, outputFile, packSize);

    return EXIT_SUCCESS;
}



