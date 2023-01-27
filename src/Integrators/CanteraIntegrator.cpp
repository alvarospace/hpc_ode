#include <vector>
#include <stdexcept>
#include <sstream>

#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Integrators/CanteraIntegrator.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "cantera/thermo.h"
#include "cantera/zerodim.h"

using std::vector;

void CanteraIntegrator::init(std::shared_ptr<Context> _ctx, IntegratorConfig config) {
    Integrator::init(_ctx, config);
    
    // Check compatibility
    auto solution = Cantera::newSolution(mechanism);
    auto thermo = solution->thermo();
    int nsp = thermo->nSpecies();
    solution.reset();
    thermo.reset();

    if (nsp != systemSize) {
        std::stringstream ss;
        ss << "mechanism size of: " << mechanism;
        ss << " (" << nsp << ")";
        ss << " differs from the mesh size: " << systemSize;
        logger->error(ss.str());
        throw std::runtime_error(ss.str());
    }
    logger->info("CanteraIntegrator initialized");
}

void CanteraIntegrator::integrate(double t0, double t) {
    Integrator::integrate(t0, t);
    double* temperatures = nullptr;
    double* enthalpies = nullptr;
    double* species = nullptr;

    if (mesh->hasEnthalpy())
        enthalpies = mesh->getEnthalpyPointer();

    temperatures = mesh->getTemperaturePointer();

    for (int i = 0; i < totalSize; i++) {
        species = mesh->getSpeciesPointer(i);
        if (mesh->hasEnthalpy())
            integrateSystem(temperatures[i], enthalpies[i], species, t);
        else
            integrateSystem(temperatures[i], species, t);
    }
}

void CanteraIntegrator::integrateSystem(double& temperature, double* species, double dt) {

    auto solution = Cantera::newSolution(mechanism);
    auto thermo = solution->thermo();
    Cantera::IdealGasConstPressureReactor reactor;
    Cantera::ReactorNet canteraIntegrator;

    thermo->setState_TPY(temperature, pressure, species);
    reactor.insert(solution);

    canteraIntegrator.setInitialTime(0.0);
    canteraIntegrator.setTolerances(reltol, abstol);
    canteraIntegrator.addReactor(reactor);

    canteraIntegrator.advance(dt);

    // Take solution
    temperature = reactor.temperature();
    thermo->getMassFractions(species);
}

void CanteraIntegrator::integrateSystem(double& temperature, double& enthalpy, double* species, double dt) {

    auto solution = Cantera::newSolution(mechanism);
    auto thermo = solution->thermo();
    Cantera::IdealGasConstPressureReactor reactor;
    Cantera::ReactorNet canteraIntegrator;

    thermo->setState_TPY(temperature, pressure, species);
    reactor.insert(solution);

    canteraIntegrator.setInitialTime(0.0);
    canteraIntegrator.setTolerances(reltol, abstol);
    canteraIntegrator.addReactor(reactor);

    canteraIntegrator.advance(dt);

    // Take solution
    temperature = reactor.temperature();
    enthalpy = thermo->enthalpy_mass();
    thermo->getMassFractions(species);
}

