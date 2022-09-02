#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"
#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

#include <string>

using namespace std;

IntegratorConfig setup(string filename, string mechanism) {
    csvReader reader(filename);
    reader.read();

    Mesh& mesh = Mesh::get();

    IntegratorConfig config;
    config.reltol = 1.0e-6;
    config.abstol = 1.0e-10;
    config.mechanism = mechanism;
    config.pressure = 101325.15;
    config.systemSize = mesh.numSpecies();
    config.totalSize = mesh.totalSize();

    return config;
}

void clean() {
    Mesh& mesh = Mesh::get();
    mesh.clear();
}

void dummy() {
    IntegratorConfig config = setup("data/res_gri_100.csv", "gri30.yaml");
    CVodeIntegrator integrator;
    integrator.init(config);
    integrator.integrate(0,1e-3);

    clean();
}

int main() {
    dummy();
    return 0;
}