#include <vector>

#include "ODEIntegrator/Integrators/CVodeIntegratorOMP.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

#include "Integrators/OpenMPRuntime.hpp"

using std::vector;

void CVodeIntegratorOMP::init(std::shared_ptr<Context> ctx, IntegratorConfig config) {
    CVodeIntegrator::init(ctx, config);
    setOMPRuntime(config.ompConfig, logger);
    logger->info("CVodeIntegratorOMP initialized");
}

void CVodeIntegratorOMP::integrate(double t0, double t) {
    Integrator::integrate(t0, t);
    // Systems allocation
    vector<vector<double>> systemsData = data_transfer_from_mesh();

    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < totalSize; i++) {
        integrateSystem(systemsData[i].data(), t);
    }

    data_transfer_to_mesh(systemsData);
}

vector<vector<double>> CVodeIntegratorOMP::data_transfer_from_mesh() {

    vector<vector<double>> systemsData(totalSize, vector<double>(systemSize,0.0f));

    vector<double> temperatures = mesh->getTemperatureVector();
    
    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < totalSize; i++) {
        vector<double> species;
        species = mesh->getSpeciesVector(i);
        systemsData[i][0] = temperatures[i];
        copy(begin(species), end(species) - 1, begin(systemsData[i]) + 1);
    }

    return systemsData;
}

void CVodeIntegratorOMP::data_transfer_to_mesh(vector<vector<double>> systemsData) {
    double* temperatures = mesh->getTemperaturePointer();
    double* species;

    #pragma omp parallel for schedule(runtime) private(species)
    for (int i = 0; i < totalSize; i++) {
        temperatures[i] = systemsData[i][0];
        species = mesh->getSpeciesPointer(i);
        copy(begin(systemsData[i]) + 1, end(systemsData[i]), species);
        species[systemSize - 1] = last_specie_calculation(species);
    }
}