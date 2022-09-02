#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

// C code from mechanism
#ifdef __cplusplus
extern "C" {
#endif
#include "Mechanism/CPU/header.h"
#ifdef __cplusplus
}
#endif

#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

int dummy_func1(double a, N_Vector b, N_Vector c, void* d) {return 0;}
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
    vector<vector<double>> systems_data = data_transfer_from_mesh(mesh);

    for (int i = 0; i < totalSize; i++) {
        
    }
}

void CVodeIntegrator::clean() {}


vector<vector<double>> CVodeIntegrator::data_transfer_from_mesh(Mesh& mesh) {

    vector<vector<double>> systems_data(totalSize, vector<double>(systemSize,0.0f));

    vector<double> temperatures = mesh.getTemperatureVector();
    vector<double> species;
    for (int i = 0; i < totalSize; i++) {
        species = mesh.getSpeciesVector(i);
        systems_data[i][0] = temperatures[i];
        move(begin(species), end(species) - 1, begin(systems_data[i]) + 1);
    }

    return move(systems_data);
}

dydt_driver CVodeIntegrator::dydt_func() {
    return dummy_func1;
}
jacobian_driver CVodeIntegrator::jacobian_func() {
    return dummy_func2;
}