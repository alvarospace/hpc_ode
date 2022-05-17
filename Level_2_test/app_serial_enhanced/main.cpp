// CVODE INCLUDES
#include <app_serial_enhanced/cvode_user.hpp>

// UTILS
#include <app_serial_enhanced/utils.hpp>

// Include for memcpy
#include <cstring>


void cvode_run(const std::string& inputFile, const std::string& outputFile,
                const size_t pack_size, const bool log) {

    /********** INPUT CONSTANTS ************/

    // Pressure Constant (Pa)
    realtype P = 101325.15;

    // Time step
    realtype t0 = 0.0f;
    realtype dt = 1;

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
        SUNLogger_SetErrorFilename(logger, "stderr.out");
        SUNLogger_SetWarningFilename(logger, "stderr.out");
        SUNLogger_SetInfoFilename(logger, "cvode_analytic_sys.info.log");
    }

    SUNDIALS_MARK_FUNCTION_BEGIN(profobj);
    
    for (int i = 0; i < n_size; i++) {
        SUNDIALS_MARK_BEGIN(profobj, "Setup");

        std::cout << "Iteration number: " << i << std::endl;

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
        retval = CVodeInit(cvode_mem, dydt_cvode, t0, y);
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

        SUNDIALS_MARK_END(profobj, "Setup");

        /***** INTEGRATION *******/
        SUNDIALS_MARK_BEGIN(profobj, "Integration");

        simTime.tic();
        retval = CVode(cvode_mem, dt, y, &t0, CV_NORMAL);
        
        // Calculation of the last element (mass convervation)
        realtype y_nsp = 0.0f;
        for (int j = 1; j < nsp; j++) {
            y_nsp += yptr[j];
        }
        y_nsp = 1 - y_nsp;

        simTime.toc();
        SUNDIALS_MARK_END(profobj, "Integration");
        /*************************/


        std::cout << "retval CVode: " << retval << std::endl;
        std::cout << "Point: " << i << "\tTime: " << simTime.time() << std::endl;

        /* Save mass fraction results */
        mesh->temp[i] = yptr[0];
        std::memcpy(mesh->matSp[i].data(), yptr+1, (nsp-1) * sizeof(realtype));
        mesh->matSp[i].data()[nsp-1] = y_nsp;

        /* Free Memory */
        N_VDestroy(y);
        CVodeFree(&cvode_mem);                   
        SUNLinSolFree(LS);                        
        SUNMatDestroy(J);
    }
    SUNDIALS_MARK_FUNCTION_END(profobj);
}

int main(int argc, char *argv[]) {
    std::string inputFile {"/home/almousa/TFM/hpc_cvode/ref_data/res_gri_1.csv"}; 
    std::string outputFile {"out.csv"};
    size_t pack_size {10};
    bool log {true};

    cvode_run(inputFile, outputFile, pack_size, log);
}



