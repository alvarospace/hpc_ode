#include <memory>
#include <algorithm>
#include <stdexcept>

#include "ODEIntegrator/Integrators/MechanismDriver.hpp"
#include "ODEIntegrator/Integrators/CVodeIntegratorGPU.hpp"
#include "ODEIntegrator/Integrators/Drivers/CVodeGPUDriver.cuh"
#include "ODEIntegrator/Integrators/CVodeGPUDataModels.hpp"
#include "Mechanism/GPU/mechanism.cuh"

#include "sundials/sundials_nvector.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_context.h"
#include "sundials/sundials_linearsolver.h"
#include "cvode/cvode.h"

#include "nvector/nvector_cuda.h"
#include "sunmatrix/sunmatrix_magmadense.h"
#include "sunlinsol/sunlinsol_magmadense.h"

void CVodeIntegratorGPU::init(std::shared_ptr<Context> ctx, IntegratorConfig config) {
    // Do not call CVodeIntegrator init method but Integrator init
    Integrator::init(ctx, config);
    uData = std::make_unique<GPUUserData>();
    uData->pressure = pressure;
    uData->ctx = ctx;
    uData->pyjac_mem = std::make_unique<mechanism_memory>();
}

void CVodeIntegratorGPU::integrate(double t0, double t) {
    Integrator::integrate(t0, t);

    int processed_systems = 0;
    while (processed_systems < totalSize) {
        int systems_to_process;
        int padded = calc_gpu_points(ctx, totalSize - processed_systems, systems_to_process);

        integrateBatch(t, processed_systems, systems_to_process, padded);

        processed_systems += systems_to_process;
    }
}

void CVodeIntegratorGPU::integrateBatch(double dt, int processed_systems, int systems_to_process, int padded) {
    // Sundials/cvode general objects
    sundials::Context sunctx;
    N_Vector y;
    SUNMatrix J;
    SUNLinearSolver LS;
    SUNMemoryHelper cuda_mem_help;
    cuda_mem_help = SUNMemoryHelper_Cuda(sunctx);
    void* cvode_mem;
    int maxSteps {10000};
    double t0 = 0.0f;
    int retVal;

    init_memory_gpu(ctx, padded, uData->pyjac_mem.get());

    // Problem size for GPU
    sunindextype nsp_GPU = NSP * systems_to_process;
    uData->nEquations = nsp_GPU;
    uData->nSystems = systems_to_process;
    try {
        // Vector allocation and initial conditions
        y = N_VNew_Cuda(nsp_GPU, sunctx);
        check_return_value(static_cast<void*>(y), "N_VNew_CUDA", 0);
        
        double* yptr = N_VGetHostArrayPointer_Cuda(y);
        setGPUArrayInitialConditions(yptr, processed_systems, systems_to_process);
        N_VCopyToDevice_Cuda(y);

        /* Call CVodeCreate to create the solver memory and specify the
        * Backward Differentiation Formula */
        cvode_mem = CVodeCreate(CV_BDF, sunctx);
        check_return_value(static_cast<void*>(cvode_mem), "CVodeCreate", 0);
        
        /* CVode init dydt = dydt_cvode(t0, y) */
        retVal = CVodeInit(cvode_mem, this->dydt_func(), t0, y);
        check_return_value(&retVal, "CVodeInit", 1);

        /* Call CVodeSStolerances */
        retVal = CVodeSStolerances(cvode_mem, reltol, abstol);
        check_return_value(&retVal, "CVodeSStolerances", 1);

        /* Create SUNMatrix for use in linear solves */
        J = SUNMatrix_MagmaDenseBlock(systems_to_process, NSP, NSP, SUNMEMTYPE_DEVICE, cuda_mem_help, NULL, sunctx);
        check_return_value(static_cast<void*>(J), "SUNMatrix_MagmaDenseBlock", 0);

        /* Create dense SUNLinearSolver (for magma library) */
        LS = SUNLinSol_MagmaDense(y, J, sunctx);
        check_return_value(static_cast<void*>(LS), "SUNLinSol_MagmaDense", 0);

        /* Attach the matrix and linear solver */
        retVal = CVodeSetLinearSolver(cvode_mem, LS, J);
        check_return_value(&retVal, "CVodeSetLinearSolver", 1);

        /* Set the user-supplied Jacobian routine Jac */
        retVal = CVodeSetJacFn(cvode_mem, this->jacobian_func());
        check_return_value(&retVal, "CVodeSetJacFn", 1);

        /* Max num steps of the method */
        retVal = CVodeSetMaxNumSteps(cvode_mem, maxSteps);
        check_return_value(&retVal, "CVodeSetMaxNumSteps", 1);

        retVal = CVodeSetUserData(cvode_mem, uData.get());
        check_return_value(&retVal, "CVodeSetUserData", 1);

        CVode(cvode_mem, dt, y, &t0, CV_NORMAL);

        // Recover and save results
        N_VCopyFromDevice_Cuda(y);
        yptr = N_VGetHostArrayPointer_Cuda(y);

        saveGPUArrayResults(yptr, processed_systems, systems_to_process);

    } catch (std::runtime_error const& e) {
        logger->error(e.what());
        free_memory_gpu(ctx, uData->pyjac_mem.get());
    }

    free_memory_gpu(ctx, uData->pyjac_mem.get());
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    N_VDestroy(y);
    SUNMatDestroy(J);
}

void CVodeIntegratorGPU::setGPUArrayInitialConditions(double *yptr, int processed_systems, int systems_to_process) {
    double *tempPtr = mesh->getTemperaturePointer();
    for (int j = 0; j < systems_to_process; j++) {
        // Index of the global point
        size_t gIndex = processed_systems + j;
        size_t lIndex = j * NSP;

        double *spPtr = mesh->getSpeciesPointer(gIndex);
        
        yptr[lIndex] = tempPtr[gIndex];
        std::copy(spPtr, spPtr + (NSP - 1), &yptr[lIndex] + 1);
    }
}

void CVodeIntegratorGPU::saveGPUArrayResults(double *yptr, int processed_systems, int systems_to_process) {
    double *tempPtr = mesh->getTemperaturePointer();
    for (int j = 0; j < systems_to_process; j++) {
        // Index of the global point
        size_t gIndex = processed_systems + j;
        size_t lIndex = j * NSP;

        realtype aux = 0.0f;
        for (int k = 1; k < NSP; k++) {
            aux += yptr[lIndex + k];
        }
        double result = 1.0f - aux;

        double *spPtr = mesh->getSpeciesPointer(gIndex);

        tempPtr[gIndex] = yptr[lIndex];
        std::copy(&yptr[lIndex] + 1, &yptr[lIndex] + 1 + (NSP - 1), spPtr);
        spPtr[NSP-1] = result;
    }
}

void CVodeIntegratorGPU::check_return_value(void* returnValue, std::string const funcName, int const opt) {
    int* retVal;
    std::stringstream ss;
    if (opt == 0 && returnValue == NULL) {
        ss << "SUNDIALS ERROR: " << funcName << "() failed - returned NULL pointer";
        throw std::runtime_error(ss.str());

    } else if (opt == 1) {
        retVal = static_cast<int*>(returnValue);
        if (*retVal < 0) {
            ss << "SUNDIALS ERROR: " << funcName << "() failed with returned value = " << *retVal;
            throw std::runtime_error(ss.str());
        }

    } else if (opt == 2 && returnValue == NULL) {
        ss << "MEMORY ERROR " << funcName << "() failed- - returned NULL pointer";
        throw std::runtime_error(ss.str()); 
    }
}

dydt_driver CVodeIntegratorGPU::dydt_func() {
    return dydt_cvode_GPU;
}

jacobian_driver CVodeIntegratorGPU::jacobian_func() {
    return jacobian_cvode_GPU;
}

