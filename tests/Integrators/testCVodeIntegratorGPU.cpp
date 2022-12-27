#include <string>
#include <memory>
#include <algorithm>

#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/Context/OutFileService.hpp"
#include "ODEIntegrator/Timer/Timer.hpp"
#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Integrators/CVodeIntegratorGPU.hpp"
// #include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"

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

void testSerialvsGPU() {
    string const FILENAME {"data/res_gri_100.csv"};
    string const MECHANISM {"gri30.yaml"};
    double const dt {1e-3};
    Timer serialTimer, GPUTimer;

    auto fileService = make_shared<OutFileService>();
    auto ctx = std::make_shared<Context>(fileService);
    auto mesh = ctx->getMesh();

    IntegratorConfig config = setup(FILENAME, MECHANISM, ctx);

    // Serial
    // serialTimer.tic();
    // CVodeIntegrator serialIntegrator;
    // serialIntegrator.init(ctx, config);
    // serialIntegrator.integrate(0,dt);
    // serialTimer.toc();
    // vector<double> voutSerial = mesh->getSpeciesVector(0);
    // mesh->clear();

    // GPU
    // config = setup(FILENAME, MECHANISM, ctx);

    GPUTimer.tic();
    CVodeIntegratorGPU GPUIntegrator;
    GPUIntegrator.init(ctx, config);
    GPUIntegrator.integrate(0, dt);
    GPUTimer.toc();
    vector<double> voutGPU = mesh->getSpeciesVector(0);
    mesh->clear();

    // // Compare
    // bool speciesEqual = equal(begin(voutSerial), end(voutSerial), begin(voutGPU), [](double const& i, double const& j) {
    //     double const err = 1e-70;
    //     return abs(i - j) <= err;
    // });
    // assert(speciesEqual);
}

int main() {
    testSerialvsGPU();
    return 0;
}