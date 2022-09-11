// CVODE INCLUDES
#include <app_magma/cvode_user.cuh>

// UTILS
#include <libutils/utils.hpp>

// Include for memcpy
#include <cstring>

// Testing
#include <app_magma/testing.hpp>

// For smart pointers
#include <memory>


void cvode_run(const std::string& inputFile, const std::string& outputFile) {

    /********** INPUT CONSTANTS ************/
    // Pressure Constant (Pa)
    realtype P = 101325.15;

    // Time step
    realtype dt = 1e-3;

    // Tolorences 
    realtype reltol = RCONST(1.0e-6);
    realtype abstol = RCONST(1.0e-10);
    /****************************/

    /*********** READ CSV FILE and SET PACKAGE SIZE ************/
    std::shared_ptr<Utils::ThermData> mesh = Utils::readCsv(inputFile);
    mesh->info();
    
    // Number of points of the mesh
    size_t n_size = mesh->points;
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

    // Problem size has to be "sunindextype", NSP * nBatch
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

    // TIMER
    Utils::Timer simTime;
    simTime.tic();

    size_t calculated_points = 0;
    while (calculated_points < n_size) {

        int gpu_points;
        int padded = calc_gpu_points(n_size - calculated_points, gpu_points);
        /* Start initialization rutine for GPU */

        // GPU PyJac memory allocation requirements
        mechanism_memory *pyjac_mem;
        //mechanism_memory *d_mem;
        pyjac_mem = new mechanism_memory;
        //initialize_gpu_memory(padded, &pyjac_mem, &d_mem);
        init_memory_gpu(padded, pyjac_mem);
        

        // CVode User defined data to access in user supplied functions
        UserData *userData = new UserData;
        userData->Pressure = P;
        userData->pyjac_mem = pyjac_mem;
        //userData->d_mem = d_mem;


        // Problem size for GPU
        nsp_GPU = NSP * gpu_points;
        userData->nEquations = nsp_GPU;
        userData->nSystems = gpu_points;

        // Vector allocation and initial conditions
        y = N_VNew_Cuda(nsp_GPU, sunctx);
        if (check_retval((void *)y, "N_VNew_CUDA", 0)) exit(EXIT_FAILURE);
        yptr = N_VGetHostArrayPointer_Cuda(y);

        for (int j = 0; j < gpu_points; j++) {
            // Index of the global point
            size_t gIndex = calculated_points + j;
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

        #ifdef TESTING
        /* TESTING "y", "dy" and Jacobian arrays */
        // Pyjac host
        std::unique_ptr<double> yPyjacHOST(new double[padded * NSP]);
        std::unique_ptr<double> dyPyjacHOST(new double[padded * NSP]);
        std::unique_ptr<double> JacPyjacHOST(new double[padded * NSP * NSP]);

        // Sundials GPU
        double *yptrGPU = N_VGetDeviceArrayPointer_Cuda(y);
        double *JacSunGPU = SUNMatrix_MagmaDense_Data(J);

        // Sundials HOST
        std::unique_ptr<double> JacSunHost(new double[gpu_points * NSP * NSP]);
        
        // Debugger objects
        std::unique_ptr<Testing::YsunYpyjac> test_y_sun_vs_py(new Testing::YsunYpyjac(
                        (double*) yptr, yptrGPU, yPyjacHOST.get(), pyjac_mem->y,
                        dyPyjacHOST.get(), pyjac_mem->dy));

        std::unique_ptr<Testing::JacSunJacPyjac> test_jacobian(new Testing::JacSunJacPyjac(
                        JacSunHost.get(), JacSunGPU, JacPyjacHOST.get(), pyjac_mem->jac));

        test_y_sun_vs_py->set_simulation_size(gpu_points, NSP, padded);
        test_y_sun_vs_py->set_logger_level("error");

        test_jacobian->set_simulation_size(gpu_points, NSP, padded);
        test_jacobian->set_logger_level("info");
        userData->test_y_sun_vs_py = test_y_sun_vs_py.get();
        userData->test_jacobian = test_jacobian.get();
        /**********************/
        #endif

        /****************** Run simulation *********************/
        realtype t0 = 0.0f;
        retval = CVode(cvode_mem, dt, y, &t0, CV_NORMAL);
        N_VCopyFromDevice_Cuda(y);
        /*******************************************************/

        /*
        *   - Calculation of the last element (mass convervation) of every point
        *   - Save data in mesh structure
        */
        std::vector<double> yLast(gpu_points, 0.0f);
        for (int j = 0; j < gpu_points; j++) {

            // Index of the global point
            size_t gIndex = calculated_points + j;
            size_t lIndex = j * NSP;

            realtype aux = 0.0f;
            for (int k = 1; k < NSP; k++) {
                aux += yptr[lIndex + k];
            }
            yLast[j] = 1.0f - aux;

            mesh->temp[gIndex] = yptr[lIndex];
            std::memcpy(mesh->matSp[gIndex].data(), &yptr[lIndex] + 1, (NSP-1) * sizeof(realtype));
            mesh->matSp[gIndex][NSP-1] = yLast[j];
        }


        /* Free Memory */
        //free_gpu_memory(&pyjac_mem, &d_mem);
        free_memory_gpu(pyjac_mem);
        delete(pyjac_mem);
        delete(userData);

        N_VDestroy(y);
        CVodeFree(&cvode_mem);
        SUNLinSolFree(LS);
        SUNMatDestroy(J);

        // Calculated points in this iteration
        calculated_points += gpu_points;
    }
    simTime.toc();

    std::cout << "Calculation time: " << simTime.time() << " s" << std::endl;

    /* Write results for validation */
    Utils::writeCsv(mesh, outputFile); 
}

int main(int argc, char *argv[]) {
    std::string inputFile {"../../../ref_data/res_gri_32.csv"}; 
    std::string outputFile {"results.csv"};

     /* Command-line arguments logic */
    std::vector<std::string> args;
    if (argc > 1) {
        for (int i = 0; i < argc; i++) {
            args.push_back(std::string(argv[i]));
        }

        if (args[1].compare("--help") == 0 || args[1].compare("-h") == 0) {
            std::cout << "Usage: " << args[0] <<" [input file=../../../ref_data/res_gri3.0.csv] [output file=results.csv]" << std::endl;
            return 0;
        } else {
            if (argc == 2) {
                inputFile = args[1];
            } else if (argc == 3) {
                inputFile = args[1];
                outputFile = args[2];
            } else {
                std::cout << "Too many arguments..." << std::endl;
                return EXIT_FAILURE;
            }
        }
    }

    /* Run integrator */
    cvode_run(inputFile, outputFile);
    
    return EXIT_SUCCESS;
}



