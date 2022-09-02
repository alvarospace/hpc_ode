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

// C code from mechanism
#ifdef __cplusplus
extern "C" {
#endif
#include "Mechanism/CPU/header.h"
#include "Mechanism/CPU/dydt.h"
#ifdef __cplusplus
}
#endif

#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

class userData {
    public:
        double pressure;
};


int dydt_cvode(double t, N_Vector y, N_Vector ydot, void* userdata) {
    double* yptr = N_VGetArrayPointer(y);
    double* ydotptr = N_VGetArrayPointer(ydot);
    userData* udata = static_cast<userData*>(userdata);

    dydt(t, udata->pressure, yptr, ydotptr);

    return 0;
}
int dummy_func2(double a, N_Vector b, N_Vector c, SUNMatrix d,  void* e, N_Vector f, N_Vector g, N_Vector h) {return 0;}

void CVodeIntegrator::init(IntegratorConfig config) {
    Integrator::init(config);

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

void CVodeIntegrator::clean() {}


vector<vector<double>> CVodeIntegrator::data_transfer_from_mesh(Mesh& mesh) {

    vector<vector<double>> systemsData(totalSize, vector<double>(systemSize,0.0f));

    vector<double> temperatures = mesh.getTemperatureVector();
    vector<double> species;
    for (int i = 0; i < totalSize; i++) {
        species = mesh.getSpeciesVector(i);
        systemsData[i][0] = temperatures[i];
        move(begin(species), end(species) - 1, begin(systemsData[i]) + 1);
    }

    return move(systemsData);
}

void CVodeIntegrator::data_transfer_to_mesh(Mesh& mesh, vector<vector<double>> systemsData) {
    double* temperatures = mesh.getTemperaturePointer();
    double* species;
    for (int i = 0; i < totalSize; i++) {
        temperatures[i] = systemsData[i][0];
        species = mesh.getSpeciesPointer(i);
        move(begin(systemsData[i]) + 1, end(systemsData[i]), species);
    }
}

void CVodeIntegrator::integrateSystem(double* system, double dt) {
    // Sundials/cvode general objects
    sundials::Context sunctx;
    N_Vector y;
    SUNMatrix J;
    SUNLinearSolver LS;
    void *cvode_mem;

    y = N_VNew_Serial(systemSize, sunctx);
    double* yptr = N_VGetArrayPointer(y);
    copy(system, system + systemSize, yptr);

    for (int i = 0; i < systemSize; i++) {
        cout << "yprt[" << i << "]: " << yptr[i] << " ";
        cout << "system[" << i << "]: " << system[i] << endl;
    }

    CVodeCreate(CV_BDF, sunctx);
    CVodeInit(cvode_mem, this->dydt_func(), 0.0, y);
    

    N_VDestroy(y);
}

dydt_driver CVodeIntegrator::dydt_func() {
    return dydt_cvode;
}
jacobian_driver CVodeIntegrator::jacobian_func() {
    return dummy_func2;
}