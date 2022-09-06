#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"
#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Integrators/CanteraIntegrator.hpp"
#include "ODEIntegrator/Integrators/CVodeIntegratorOMP.hpp"
#include "ODEIntegrator/Timer/Timer.hpp"

#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>

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
    string filename {"data/res_gri_1.csv"};
    string mechanism {"gri30.yaml"};

    IntegratorConfig config = setup(filename, mechanism);
    double dt = 1e-3;
    Mesh& mesh = Mesh::get();

    double tin = mesh.getTemperatureVector()[0];
    vector<double> vin = mesh.getSpeciesVector(0);
    
    //CVode
    CVodeIntegrator cvodeIntegrator;
    cvodeIntegrator.init(config);
    cvodeIntegrator.integrate(0,dt);
    double tempCVode = mesh.getTemperatureVector()[0];
    vector<double> voutCVode = mesh.getSpeciesVector(0);
    clean();

    config = setup(filename, mechanism);
    //Cantera
    CanteraIntegrator canteraIntegrator;
    cvodeIntegrator.init(config);
    cvodeIntegrator.integrate(0,dt);
    double tempCantera = mesh.getTemperatureVector()[0];
    vector<double> voutCantera = mesh.getSpeciesVector(0);
    clean();

    // Species comparison
    bool speciesEqual = equal(begin(voutCVode), end(voutCVode), begin(voutCantera), [](double const& i, double const& j) {
        double const err = 1e-70;
        return abs(i - j) <= err;
    });
    assert(speciesEqual);
}

void testSerialvsOMP() {
    string filename {"data/res_gri_100.csv"};
    string mechanism {"gri30.yaml"};
    double dt = 1e-3;
    Timer serialTimer, OMPTimer;

    Mesh& mesh = Mesh::get();

    IntegratorConfig config = setup(filename, mechanism);

    // Serial
    serialTimer.tic();
    CVodeIntegrator serialIntegrator;
    serialIntegrator.init(config);
    serialIntegrator.integrate(0,dt);
    serialTimer.toc();
    vector<double> voutSerial = mesh.getSpeciesVector(0);
    clean();

    config = setup(filename, mechanism);

    // OMP
    OMPTimer.tic();
    CVodeIntegratorOMP OMPIntegrator;
    OMPIntegrator.init(config);
    OMPIntegrator.integrate(0, dt);
    OMPTimer.toc();
    vector<double> voutOMP = mesh.getSpeciesVector(0);
    clean();

    bool speciesEqual = equal(begin(voutSerial), end(voutSerial), begin(voutOMP), [](double const& i, double const& j) {
        double const err = 1e-70;
        return abs(i - j) <= err;
    });
    assert(speciesEqual);

    // OMP takes less time than serial integration
    assert(serialTimer.getTime() > OMPTimer.getTime());
    cout << serialTimer.getTime() << endl;
    cout << OMPTimer.getTime() << endl;
}

int main() {
    testCVodeIntegrator();
    testCVodevsCantera();
    testSerialvsOMP();
    return 0;
}