// CVODE INCLUDES
#include <app_serial_enhanced/cvode_user.hpp>

// UTILS
#include <app_serial_enhanced/utils.hpp>

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
    size_t mod = n_size % pack_size;
    size_t ext_it = n_size / pack_size;
    // If the division is not exact, add 1 extra iteration (the remainder)
    if (mod > 0)
        ext_it++;

    std::cout << "Package size: " << pack_size << std::endl;
    std::cout << "Number of packages: " << ext_it << std::endl;
    std::cout << "Remainder points (will be process by the last package): " << mod << std::endl;
    /************************************************************/

    // Check that number of species read is the same that in pyjac library
    if (mesh->nsp != NSP) {
        std::cout << "Error: Number of species in the input file does not match Pyjac species" << std::endl;
        exit(EXIT_FAILURE);
    }

    /********** CVODE INTEGRATION ************/

    // Variables declarations
    int m_maxsteps = 10000;

    // Problem size has to be "sunindextype"
    sunindextype nsp = NSP;

    /***** MULTI-THREADING objects management *****/
    int num_threads = Utils::printNumThreads();

    /* Sundials context creation and automatic free memory (C++ feature) */
    SUNContext sunctxs[num_threads];

    /* Flags to verificate initialization */
    int cvode_initialized[num_threads];

    /* Cvode mem */
    void* cvode_mem[num_threads];

    /* Data variables */
    N_Vector y[num_threads];
    realtype *yptr[num_threads];
    SUNMatrix J[num_threads];
    SUNLinearSolver LS[num_threads];

    int retval[num_threads];

    /**********************************************/

    // Cvode creation
    for (int i = 0; i < num_threads; i++) {
        retval[i] = SUNContext_Create(NULL, &sunctxs[i]);
        cvode_mem[i] = CVodeCreate(CV_BDF, sunctxs[i]);
        cvode_initialized[i] = 0;
    }

    
    // Integration
    simTime.tic();

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < ext_it; i++) {
        size_t thread_points;

        // Limit the pack_size to the last iteration
        // because the number of points to process is the remainder
        if (i == (ext_it - 1) && mod > 0)
            thread_points = mod;
        else
            thread_points = pack_size;

        // Thread ID
        int id = omp_get_thread_num();

        for (int j = 0; j < thread_points; j++) {
            // Index of the global point
            size_t index = i * pack_size + j;

            // If first time integration set cvode
            if (!cvode_initialized[id]) {

                cvode_initialized[id] = 1;

                /* Memory allocation for initial conditions */
                y[id] = N_VNew_Serial(nsp, sunctxs[id]);
                if (check_retval((void *)y[id], "N_VNew_Serial", 0)) exit(EXIT_FAILURE);
                yptr[id] = NV_DATA_S(y[id]);

                /* Initials conditions of the current system */
                // Y = {Temp, Y1, Y2, ..., Y_nsp - 1} Last species can be solved with mass conservation:
                // Y_nsp = 1 - SUM_i ( Y_i )
                yptr[id][0] = mesh->temp[index];
                std::memcpy(yptr[id]+1, mesh->matSp[index].data(), (nsp-1) * sizeof(realtype));

                /* CVode init dydt = dydt_cvode(t0, y) */
                retval[id] = CVodeInit(cvode_mem[id], dydt_cvode, 0.0, y[id]);
                if (check_retval(&retval[id], "CVodeInit", 1)) exit(EXIT_FAILURE);

                /* Call CVodeSStolerances */
                retval[id] = CVodeSStolerances(cvode_mem[id], reltol, abstol);
                if (check_retval(&retval[id], "CVodeSVtolerances", 1)) exit(EXIT_FAILURE);

                /* Create dense SUNMatrix for use in Jacobian evaluation and linear solver */
                J[id] = SUNDenseMatrix(nsp, nsp, sunctxs[id]);
                if (check_retval((void *)J[id], "SUNDenseMatrix", 0)) exit(EXIT_FAILURE);

                /* Create dense SUNLinearSolver (LU decomposition) */
                LS[id] = SUNLinSol_Dense(y[id], J[id], sunctxs[id]);
                if (check_retval((void *)LS[id], "SUNLinSol_Dense", 0)) exit(EXIT_FAILURE);

                /* Attach the matrix and linear solver */
                retval[id] = CVodeSetLinearSolver(cvode_mem[id], LS[id], J[id]);
                if (check_retval(&retval[id], "CVodeSetLinearSolver", 1)) exit(EXIT_FAILURE);

                /* Set the user-supplied Jacobian routine Jac */
                retval[id] = CVodeSetJacFn(cvode_mem[id], eval_jacob_cvode);
                if (check_retval(&retval[id], "CVodeSetJacFn", 1)) exit(EXIT_FAILURE);

                /* Max num steps of the method */
                retval[id] = CVodeSetMaxNumSteps(cvode_mem[id], m_maxsteps);
                if (check_retval(&retval[id], "CVodeSetMaxNumSteps", 1)) exit(EXIT_FAILURE);

                /* Save Constant Pressure internally in CVode memory (user data for internal functions) */
                retval[id] = CVodeSetUserData(cvode_mem[id], &P);

            } else { // Reinit

                /* Initials conditions of the current system */
                yptr[id][0] = mesh->temp[index];
                std::memcpy(yptr[id]+1, mesh->matSp[index].data(), (nsp-1) * sizeof(realtype));

                retval[id] = CVodeReInit(cvode_mem[id], 0.0, y[id]);
                if (check_retval(&retval[id], "CVodeReInit", 1)) exit(EXIT_FAILURE);
            }

            // Run simulation
            realtype t0 = 0.0f;
            retval[id] = CVode(cvode_mem[id], dt, y[id], &t0, CV_NORMAL);

            // Calculation of the last element (mass convervation)
            realtype y_nsp = 0.0f;
            for (int k = 1; k < nsp; k++) {
                y_nsp += yptr[id][k];
            }
            y_nsp = 1.0 - y_nsp;

            // Retrieve results
            mesh->temp[index] = yptr[id][0];
            std::memcpy(mesh->matSp[index].data(), yptr[id]+1, (nsp-1) * sizeof(realtype));
            mesh->matSp[index][nsp-1] = y_nsp;
        }
    }

    /* Free Memory */
    for (int i = 0; i < num_threads; i++) {
        N_VDestroy(y[i]);
        CVodeFree(&cvode_mem[i]);                   
        SUNLinSolFree(LS[i]);                        
        SUNMatDestroy(J[i]);
        SUNContext_Free(&sunctxs[i]);
    }
    

    simTime.toc();

    std::cout << "Calculation time: " << simTime.time() << " s" << std::endl;

    /* Write results for validation */
    Utils::writeCsv(mesh, outputFile); 
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



