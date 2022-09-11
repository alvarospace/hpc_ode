// CVODE INCLUDES
#include <app_serial_enhanced/cvode_user.hpp>

// UTILS
#include <app_serial_enhanced/utils.hpp>

// Include for memcpy
#include <cstring>


void cvode_run(const std::string& inputFile, const std::string& outputFile, const bool log) {

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
    /************************************************************/

    // Check that number of species read is the same that in pyjac library
    if (mesh->nsp != NSP) {
        std::cout << "Error: Number of species in the input file does not match Pyjac species" << std::endl;
        exit(EXIT_FAILURE);
    }

    /********** CVODE INTEGRATION ************/

    // Variables declarations
    N_Vector y;
    realtype *yptr;
    SUNMatrix J;
    SUNLinearSolver LS;
    void *cvode_mem;
    int m_maxsteps = 10000;
    int retval;
    SUNLogger logger = nullptr;
    SUNProfiler profobj;

    // Problem size has to be "sunindextype"
    sunindextype nsp = NSP;


    /* Sundials context creation and automatic free memory (C++ feature) */
    sundials::Context sunctx = sundials::Context();
    SUNContext_GetProfiler(sunctx, &profobj);

    // LOGGER
    if (log) {
        SUNContext_GetLogger(sunctx, &logger);
        SUNLogger_SetErrorFilename(logger, "stderr");
        SUNLogger_SetWarningFilename(logger, "stderr");
        SUNLogger_SetInfoFilename(logger, "cvode_analytic_sys.info.log");
    }

    SUNDIALS_MARK_FUNCTION_BEGIN(profobj);
    simTime.tic();
    for (int i = 0; i < n_size; i++) {
        //SUNDIALS_MARK_BEGIN(profobj, "Setup");

        /* Memory allocation for initial conditions */
        y = N_VNew_Serial(nsp, sunctx);
        if (check_retval((void *)y, "N_VNew_Serial", 0)) exit(EXIT_FAILURE);
        yptr = NV_DATA_S(y);

        /* Initials conditions of the current system */
        // Y = {Temp, Y1, Y2, ..., Y_nsp - 1} Last species can be solved with mass conservation:
        // Y_nsp = 1 - SUM_i ( Y_i )
        yptr[0] = mesh->temp[i];
        std::memcpy(yptr+1, mesh->matSp[i].data(), (nsp-1) * sizeof(realtype));

        /* CVode create with Backward Differentation Formula */
        cvode_mem = CVodeCreate(CV_BDF, sunctx);
        if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) exit(EXIT_FAILURE);

        /* CVode init dydt = dydt_cvode(t0, y) */
        retval = CVodeInit(cvode_mem, dydt_cvode, 0.0, y);
        if (check_retval(&retval, "CVodeInit", 1)) exit(EXIT_FAILURE);

        /* Call CVodeSVtolerances */
        retval = CVodeSStolerances(cvode_mem, reltol, abstol);
        if (check_retval(&retval, "CVodeSVtolerances", 1)) exit(EXIT_FAILURE);

        /* Create dense SUNMatrix for use in Jacobian evaluation and linear solver */
        J = SUNDenseMatrix(nsp, nsp, sunctx);
        if (check_retval((void *)J, "SUNDenseMatrix", 0)) exit(EXIT_FAILURE);

        /* Create dense SUNLinearSolver (LU decomposition) */
        LS = SUNLinSol_Dense(y, J, sunctx);
        if (check_retval((void *)LS, "SUNLinSol_Dense", 0)) exit(EXIT_FAILURE);

        /* Attach the matrix and linear solver */
        retval = CVodeSetLinearSolver(cvode_mem, LS, J);
        if (check_retval(&retval, "CVodeSetLinearSolver", 1)) exit(EXIT_FAILURE);

        /* Set the user-supplied Jacobian routine Jac */
        retval = CVodeSetJacFn(cvode_mem, eval_jacob_cvode);
        if (check_retval(&retval, "CVodeSetJacFn", 1)) exit(EXIT_FAILURE);

        /* Max num steps of the method */
        retval = CVodeSetMaxNumSteps(cvode_mem, m_maxsteps);
        if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) exit(EXIT_FAILURE);

        /* Save Constant Pressure internally in CVode memory (user data for internal functions) */
        retval = CVodeSetUserData(cvode_mem, &P);

        //SUNDIALS_MARK_END(profobj, "Setup");

        /***** INTEGRATION *******/
        //SUNDIALS_MARK_BEGIN(profobj, "Integration");
        realtype t0 = 0.0f;
        retval = CVode(cvode_mem, dt, y, &t0, CV_NORMAL);

        // Calculation of the last element (mass convervation)
        realtype y_nsp = 0.0f;
        for (int j = 1; j < nsp; j++) {
            y_nsp += yptr[j];
        }
        y_nsp = 1.0 - y_nsp;

        //SUNDIALS_MARK_END(profobj, "Integration");
        /*************************/


        /* Save mass fraction results */
        mesh->temp[i] = yptr[0];
        std::memcpy(mesh->matSp[i].data(), yptr+1, (nsp-1) * sizeof(realtype));
        mesh->matSp[i][nsp-1] = y_nsp;

        /* Free Memory */
        N_VDestroy(y);
        CVodeFree(&cvode_mem);                   
        SUNLinSolFree(LS);                        
        SUNMatDestroy(J);
    }
    simTime.toc();
    SUNDIALS_MARK_FUNCTION_END(profobj);

    std::cout << "Calculation time: " << simTime.time() << " s" << std::endl;

    /* Write results for validation */
    Utils::writeCsv(mesh, outputFile); 

    FILE *fid;
    fid = fopen("profiling.txt", "w");
    SUNProfiler_Print(profobj, fid);
    fclose(fid);
}

int main(int argc, char *argv[]) {
    std::string inputFile {"../../../ref_data/res_gri3.0.csv"}; 
    std::string outputFile {"results.csv"};
    bool log {false};

     /* Command-line arguments logic */
    std::vector<std::string> args;
    if (argc > 1) {
        for (int i = 0; i < argc; i++) {
            args.push_back(std::string(argv[i]));
        }

        if (args[1].compare("--help") == 0 || args[1].compare("-h") == 0) {
            std::cout << "Usage: " << args[0] <<" [input file=../../../ref_data/res_gri3.0.csv] [output file=results.csv] [log=0]" << std::endl;
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
                if (std::stoi(args[3]) == 1)
                    log = true;
            } else {
                std::cout << "Too many arguments..." << std::endl;
                return EXIT_FAILURE;
            }
        }
    }

    /* Run integrator */
    cvode_run(inputFile, outputFile, log);

    return EXIT_SUCCESS;
}



