#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <memory>

#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/Context/OutFileService.hpp"
#include "ODEIntegrator/Integrators/CanteraIntegrator.hpp"
#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
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

void testCanteraIntegrator(string filename, string mechanism) {
    auto fileService = std::make_shared<OutFileService>();
    auto mesh = std::make_shared<Mesh>();
    auto ctx = std::make_shared<Context>(fileService, mesh);
    IntegratorConfig config = setup(filename, mechanism, ctx);

    vector<double> vin = mesh->getSpeciesVector(0);

    CanteraIntegrator integrator;
    integrator.init(ctx, config);
    integrator.integrate(0, 1e-3);
    
    vector<double> vout = mesh->getSpeciesVector(0);

    // Species comparison
    bool speciesEqual = equal(begin(vin), end(vin), begin(vout), [](double const& i, double const& j) {
        double const err = 1e-70;
        return abs(i - j) <= err;
    });
    assert(!speciesEqual);

    mesh->clear();
}

void testSerialvsOMP(string filename, string mechanism) {
    Timer serialTimer, OMPTimer;
    auto fileService = std::make_shared<OutFileService>();
    auto mesh = std::make_shared<Mesh>();
    auto ctx = std::make_shared<Context>(fileService, mesh);
    double dt = 1e-3;
    
    // Serial
    IntegratorConfig config = setup(filename, mechanism, ctx);
    serialTimer.tic();

    CanteraIntegrator serialIntegrator;
    serialIntegrator.init(ctx, config);
    serialIntegrator.integrate(0, dt);

    serialTimer.toc();
    vector<double> vSerial = mesh->getSpeciesVector(50);
    mesh->clear();

    // OMP
    config = setup(filename, mechanism, ctx);
    OMPTimer.tic();

    CanteraIntegratorOMP OMPIntegrator;
    OMPIntegrator.init(ctx, config);
    OMPIntegrator.integrate(0, dt);

    OMPTimer.toc();
    vector<double> vOMP = mesh->getSpeciesVector(50);

    assert(serialTimer.getTime() > OMPTimer.getTime());
    cout << serialTimer.getTime() << endl;
    cout << OMPTimer.getTime() << endl;

    // Species comparison
    bool speciesEqual = equal(begin(vSerial), end(vSerial), begin(vOMP), [](double const& i, double const& j) {
        double const err = 1e-70;
        return abs(i - j) <= err;
    });
    assert(speciesEqual);

    mesh->clear();
} 

int main() {
    testCanteraIntegrator("data/h2o2.csv", "h2o2.yaml");
    testCanteraIntegrator("data/res_gri_100.csv", "gri30.yaml");
    testSerialvsOMP("data/res_gri_100.csv", "gri30.yaml");
    return 0;
}