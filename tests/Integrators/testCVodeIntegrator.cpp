#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"
#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Integrators/CanteraIntegrator.hpp"

#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

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

void testCVodeIntegrator() {
    IntegratorConfig config = setup("data/res_gri_1.csv", "gri30.yaml");
    Mesh& mesh = Mesh::get();

    vector<double> vin = mesh.getSpeciesVector(0);

    CVodeIntegrator integrator;
    integrator.init(config);
    integrator.integrate(0,1e-3);

    vector<double> vout = mesh.getSpeciesVector(0);

    // Species comparison
    bool speciesEqual = equal(begin(vin), end(vin), begin(vout), [](double const& i, double const& j) {
        double const err = 1e-70;
        return abs(i - j) <= err;
    });
    assert(!speciesEqual);

    clean();
}

void testCVodevsCantera() {
    IntegratorConfig config = setup("data/res_gri_1.csv", "gri30.yaml");
    double dt = 1e-3;
    Mesh& mesh = Mesh::get();

    //CVode
    CVodeIntegrator cvodeIntegrator;
    cvodeIntegrator.init(config);
    cvodeIntegrator.integrate(0,dt);
    double tempCVode = mesh.getTemperatureVector()[0];
    vector<double> voutCVode = mesh.getSpeciesVector(0);
    clean();

    config = setup("data/res_gri_1.csv", "gri30.yaml");
    //Cantera
    CanteraIntegrator canteraIntegrator;
    cvodeIntegrator.init(config);
    cvodeIntegrator.integrate(0,dt);
    double tempCantera = mesh.getTemperatureVector()[0];
    vector<double> voutCantera = mesh.getSpeciesVector(0);
    clean();

    cout << "CVode\tCantera" << endl;
    cout << tempCVode << "\t" << tempCantera << endl;
    for (int i = 0; i < config.systemSize; i++) {
        cout << voutCVode[i] << "\t" << voutCantera[i] << endl;
    }
}

int main() {
    testCVodeIntegrator();
    testCVodevsCantera();
    return 0;
}