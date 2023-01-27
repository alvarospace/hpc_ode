#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>

#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"
#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Integrators/CanteraIntegrator.hpp"
#include "ODEIntegrator/Integrators/CVodeIntegratorOMP.hpp"
#include "ODEIntegrator/Timer/Timer.hpp"

using namespace std;

IntegratorConfig setup(string filename, string mechanism, shared_ptr<Context> ctx) {
    csvReader reader(ctx, filename);
    reader.read();

    IntegratorConfig config;
    config.reltol = 1.0e-6;
    config.abstol = 1.0e-10;
    config.mechanism = mechanism;
    config.pressure = 101325.15;

    return config;
}

void testCVodeIntegrator() {
    auto fileService = std::make_shared<OutFileService>();
    auto mesh = std::make_shared<Mesh>();
    auto ctx = std::make_shared<Context>(fileService, mesh);
    IntegratorConfig config = setup("data/res_gri_1.csv", "gri30.yaml", ctx);

    vector<double> vin = mesh->getSpeciesVector(0);

    CVodeIntegrator integrator;
    integrator.init(ctx, config);
    integrator.integrate(0,1e-3);

    vector<double> vout = mesh->getSpeciesVector(0);

    // Species comparison
    bool speciesEqual = equal(begin(vin), end(vin), begin(vout), [](double const& i, double const& j) {
        double const err = 1e-70;
        return abs(i - j) <= err;
    });
    assert(!speciesEqual);

    mesh->clear();
}

void testCVodevsCantera() {
    string filename {"data/res_gri_1.csv"};
    string mechanism {"gri30.yaml"};

    auto fileService = std::make_shared<OutFileService>();
    auto mesh = std::make_shared<Mesh>();
    auto ctx = std::make_shared<Context>(fileService, mesh);

    IntegratorConfig config = setup(filename, mechanism, ctx);
    double dt = 1e-3;

    double tin = mesh->getTemperatureVector()[0];
    vector<double> vin = mesh->getSpeciesVector(0);
    
    //CVode
    CVodeIntegrator cvodeIntegrator;
    cvodeIntegrator.init(ctx, config);
    cvodeIntegrator.integrate(0,dt);
    double tempCVode = mesh->getTemperatureVector()[0];
    vector<double> voutCVode = mesh->getSpeciesVector(0);
    mesh->clear();

    config = setup(filename, mechanism, ctx);
    //Cantera
    CanteraIntegrator canteraIntegrator;
    cvodeIntegrator.init(ctx, config);
    cvodeIntegrator.integrate(0,dt);
    double tempCantera = mesh->getTemperatureVector()[0];
    vector<double> voutCantera = mesh->getSpeciesVector(0);
    mesh->clear();

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

    auto fileService = std::make_shared<OutFileService>();
    auto mesh = std::make_shared<Mesh>();
    auto ctx = std::make_shared<Context>(fileService, mesh);

    IntegratorConfig config = setup(filename, mechanism, ctx);

    // Serial
    serialTimer.tic();
    CVodeIntegrator serialIntegrator;
    serialIntegrator.init(ctx, config);
    serialIntegrator.integrate(0,dt);
    serialTimer.toc();
    vector<double> voutSerial = mesh->getSpeciesVector(0);
    mesh->clear();

    config = setup(filename, mechanism, ctx);

    // OMP
    OMPTimer.tic();
    CVodeIntegratorOMP OMPIntegrator;
    OMPIntegrator.init(ctx, config);
    OMPIntegrator.integrate(0, dt);
    OMPTimer.toc();
    vector<double> voutOMP = mesh->getSpeciesVector(0);
    mesh->clear();

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