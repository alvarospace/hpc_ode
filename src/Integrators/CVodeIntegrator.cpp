#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

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
#include "Mechanism/CPU/dydt.h"
#include "Mechanism/CPU/jacob.h"
#ifdef __cplusplus
}
#endif

#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

using namespace std;

int dydt_cvode_serial(double t, N_Vector y, N_Vector ydot, void* userdata) {
    double* yptr = N_VGetArrayPointer(y);
    double* ydotptr = N_VGetArrayPointer(ydot);
    CVodeIntegrator::UserData* udata = static_cast<CVodeIntegrator::UserData*>(userdata);

    dydt(t, udata->pressure, yptr, ydotptr);

    return 0;
}
int jacobian_cvode_serial(double t, N_Vector y, N_Vector ydot, SUNMatrix J, void* userdata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    double* yptr = N_VGetArrayPointer(y);
    double* Jptr = SM_DATA_D(J);
    CVodeIntegrator::UserData* udata = static_cast<CVodeIntegrator::UserData*>(userdata);

    eval_jacob(t, udata->pressure, yptr, Jptr);

    return 0;
}

void CVodeIntegrator::init(IntegratorConfig config) {
    Integrator::init(config);
    uData = make_unique<UserData>();
    uData->pressure = pressure;

    // Check compatibility with the mechanism
    if (systemSize != NSP) {
        std::stringstream ss;
        ss << "mechanism size of: " << mechanism;
        ss << " (" << NSP << ")";
        ss << " differs from the mesh size: " << systemSize;
        throw std::logic_error(ss.str());
    }
}

void CVodeIntegrator::integrate(double t0, double t) {
    Mesh& mesh = Mesh::get();

    // Systems allocation
    vector<vector<double>> systemsData = data_transfer_from_mesh(mesh);

    for (int i = 0; i < totalSize; i++) {
        integrateSystem(systemsData[i].data(), t);
    }

    data_transfer_to_mesh(mesh, systemsData);
}

vector<vector<double>> CVodeIntegrator::data_transfer_from_mesh(Mesh& mesh) {

    vector<vector<double>> systemsData(totalSize, vector<double>(systemSize,0.0f));

    vector<double> temperatures = mesh.getTemperatureVector();
    vector<double> species;
    for (int i = 0; i < totalSize; i++) {
        species = mesh.getSpeciesVector(i);
        systemsData[i][0] = temperatures[i];
        copy(begin(species), end(species) - 1, begin(systemsData[i]) + 1);
    }

    return move(systemsData);
}

void CVodeIntegrator::data_transfer_to_mesh(Mesh& mesh, vector<vector<double>> systemsData) {
    double* temperatures = mesh.getTemperaturePointer();
    double* species;
    for (int i = 0; i < totalSize; i++) {
        temperatures[i] = systemsData[i][0];
        species = mesh.getSpeciesPointer(i);
        vector<double> species2 = mesh.getSpeciesVector(i);
        double tam = species2.size();
        move(begin(systemsData[i]) + 1, end(systemsData[i]), species);
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

    y = N_VNew_Serial(systemSize, sunctx);
    double* yptr = N_VGetArrayPointer(y);
    copy(system, system + systemSize, yptr);

    cvode_mem = CVodeCreate(CV_BDF, sunctx);
    CVodeInit(cvode_mem, this->dydt_func(), t0, y);
    CVodeSStolerances(cvode_mem, reltol, abstol);
    J = SUNDenseMatrix(systemSize, systemSize, sunctx);
    LS = SUNLinSol_Dense(y, J, sunctx);
    CVodeSetLinearSolver(cvode_mem, LS, J);
    CVodeSetJacFn(cvode_mem, this->jacobian_func());
    CVodeSetUserData(cvode_mem, uData.get());

    CVode(cvode_mem, dt, y, &t0, CV_NORMAL);

    copy(yptr, yptr + systemSize, system);

    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    N_VDestroy(y);
    SUNMatDestroy(J);
}

dydt_driver CVodeIntegrator::dydt_func() {
    return dydt_cvode_serial;
}
jacobian_driver CVodeIntegrator::jacobian_func() {
    return jacobian_cvode_serial;
}