#include "ODEIntegrator/Integrators/CanteraIntegrator.hpp"
#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Timer/Timer.hpp"

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>

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

void testCanteraIntegrator(string filename, string mechanism) {
    IntegratorConfig config = setup(filename, mechanism);

    Mesh& mesh = Mesh::get();
    vector<double> vin = mesh.getSpeciesVector(0);

    CanteraIntegrator integrator;
    integrator.init(config);
    integrator.integrate(0, 1e-3);
    
    vector<double> vout = mesh.getSpeciesVector(0);

    // Species comparison
    bool speciesEqual = equal(begin(vin), end(vin), begin(vout), [](double const& i, double const& j) {
        double const err = 1e-70;
        return abs(i - j) <= err;
    });
    assert(!speciesEqual);

    clean();
}

void testSerialvsOMP(string filename, string mechanism) {
    Timer serialTimer, OMPTimer;
    Mesh& mesh = Mesh::get();
    double dt = 1e-3;
    
    // Serial
    IntegratorConfig config = setup(filename, mechanism);
    serialTimer.tic();

    CanteraIntegrator serialIntegrator;
    serialIntegrator.init(config);
    serialIntegrator.integrate(0, dt);

    serialTimer.toc();
    vector<double> vSerial = mesh.getSpeciesVector(50);
    clean();

    // OMP
    config = setup(filename, mechanism);
    OMPTimer.tic();

    CanteraIntegratorOMP OMPIntegrator;
    OMPIntegrator.init(config);
    OMPIntegrator.integrate(0, dt);

    OMPTimer.toc();
    vector<double> vOMP = mesh.getSpeciesVector(50);

    assert(serialTimer.getTime() > OMPTimer.getTime());
    cout << serialTimer.getTime() << endl;
    cout << OMPTimer.getTime() << endl;

    // Species comparison
    bool speciesEqual = equal(begin(vSerial), end(vSerial), begin(vOMP), [](double const& i, double const& j) {
        double const err = 1e-70;
        return abs(i - j) <= err;
    });
    assert(speciesEqual);

    clean();
} 

int main() {
    testCanteraIntegrator("data/h2o2.csv", "h2o2.yaml");
    testCanteraIntegrator("data/res_gri_100.csv", "gri30.yaml");
    testSerialvsOMP("data/res_gri_100.csv", "gri30.yaml");
    return 0;
}