#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Integrators/Drivers/CVodeSerialDriver.hpp"

#include "sundials/sundials_nvector.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_context.h"
#include "sundials/sundials_linearsolver.h"
#include "cvode/cvode.h"

// Sundials serial headers
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_dense.h"

// C code from mechanism
#ifdef __cplusplus
extern "C" {
#endif
#include "Mechanism/CPU/header.h"
#ifdef __cplusplus
}
#endif

using namespace std;

void CVodeIntegrator::init(std::shared_ptr<Context> _ctx, IntegratorConfig config) {
    Integrator::init(_ctx, config);
    uData = make_unique<UserData>();
    uData->pressure = pressure;

    // Check compatibility with the mechanism
    if (systemSize != NSP) {
        std::stringstream ss;
        ss << "mechanism size of: " << mechanism;
        ss << " (" << NSP << ")";
        ss << " differs from the mesh size: " << systemSize;
        logger->error(ss.str());
        throw std::runtime_error(ss.str());
    }
    logger->info("CVodeIntegrator initialized");
}

void CVodeIntegrator::integrate(double t0, double t) {

    // Systems allocation
    vector<vector<double>> systemsData = data_transfer_from_mesh();

    for (int i = 0; i < totalSize; i++) {
        integrateSystem(systemsData[i].data(), t);
    }

    data_transfer_to_mesh(systemsData);
}

vector<vector<double>> CVodeIntegrator::data_transfer_from_mesh() {

    vector<vector<double>> systemsData(totalSize, vector<double>(systemSize,0.0f));

    vector<double> temperatures = mesh->getTemperatureVector();
    vector<double> species;
    for (int i = 0; i < totalSize; i++) {
        species = mesh->getSpeciesVector(i);
        systemsData[i][0] = temperatures[i];
        copy(begin(species), end(species) - 1, begin(systemsData[i]) + 1);
    }

    return systemsData;
}

void CVodeIntegrator::data_transfer_to_mesh(vector<vector<double>> systemsData) {
    double* temperatures = mesh->getTemperaturePointer();
    double* species;
    for (int i = 0; i < totalSize; i++) {
        temperatures[i] = systemsData[i][0];
        species = mesh->getSpeciesPointer(i);
        copy(begin(systemsData[i]) + 1, end(systemsData[i]), species);
        species[systemSize - 1] = last_specie_calculation(species);
    }
}

double CVodeIntegrator::last_specie_calculation(double* species) {
    double y_nsp {0.0f};
    for (int i = 0; i < systemSize - 1; i++) {
        y_nsp += species[i];
    }
    return 1.0f - y_nsp;
}

void CVodeIntegrator::integrateSystem(double* system, double dt) {
    // Sundials/cvode general objects
    sundials::Context sunctx;
    N_Vector y;
    SUNMatrix J;
    SUNLinearSolver LS;
    void *cvode_mem;
    double t0 = 0.0f;
    int retVal;
    int maxSteps {10000};

    try {
        y = N_VNew_Serial(systemSize, sunctx);
        check_return_value(static_cast<void*>(y), "N_VNew_Serial", 0);

        double* yptr = N_VGetArrayPointer(y);
        copy(system, system + systemSize, yptr);

        cvode_mem = CVodeCreate(CV_BDF, sunctx);
        check_return_value(static_cast<void*>(cvode_mem), "CVodeCreate", 0);

        retVal = CVodeInit(cvode_mem, this->dydt_func(), t0, y);
        check_return_value(&retVal, "CVodeInit", 1);

        retVal = CVodeSStolerances(cvode_mem, reltol, abstol);
        check_return_value(&retVal, "CVodeSStolerances", 1);

        J = SUNDenseMatrix(systemSize, systemSize, sunctx);
        check_return_value(static_cast<void*>(J), "SUNDenseMatrix", 0);

        LS = SUNLinSol_Dense(y, J, sunctx);
        check_return_value(static_cast<void*>(LS), "SUNLinSol_Dense", 0);

        retVal = CVodeSetLinearSolver(cvode_mem, LS, J);
        check_return_value(&retVal, "CVodeSetLinearSolver", 1);

        retVal = CVodeSetJacFn(cvode_mem, this->jacobian_func());
        check_return_value(&retVal, "CVodeSetJacFn", 1);

        retVal = CVodeSetMaxNumSteps(cvode_mem, maxSteps);
        check_return_value(&retVal, "CVodeSetMaxNumSteps", 1);

        retVal = CVodeSetUserData(cvode_mem, uData.get());
        check_return_value(&retVal, "CVodeSetUserData", 1);

        CVode(cvode_mem, dt, y, &t0, CV_NORMAL);

        copy(yptr, yptr + systemSize, system);
        
    } catch (runtime_error const& e) {
        logger->error(e.what());
    }

    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    N_VDestroy(y);
    SUNMatDestroy(J);
}

void CVodeIntegrator::check_return_value(void* returnValue, string const funcName, int const opt) {
    int* retVal;
    stringstream ss;
    if (opt == 0 && returnValue == NULL) {
        ss << "SUNDIALS ERROR: " << funcName << "() failed - returned NULL pointer";
        throw runtime_error(ss.str());

    } else if (opt == 1) {
        retVal = static_cast<int*>(returnValue);
        if (*retVal < 0) {
            ss << "SUNDIALS ERROR: " << funcName << "() failed with returned value = " << *retVal;
            throw runtime_error(ss.str());
        }

    } else if (opt == 2 && returnValue == NULL) {
        ss << "MEMORY ERROR " << funcName << "() failed- - returned NULL pointer";
        throw runtime_error(ss.str()); 
    }
}

dydt_driver CVodeIntegrator::dydt_func() {
    return dydt_cvode_serial;
}

jacobian_driver CVodeIntegrator::jacobian_func() {
    return jacobian_cvode_serial;
}