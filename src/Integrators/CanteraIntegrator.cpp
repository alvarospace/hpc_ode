#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Integrators/CanteraIntegrator.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "cantera/thermo.h"
#include "cantera/zerodim.h"


void CanteraIntegrator::init(IntegratorConfig config) {
    totalSize = config.totalSize;
    systemSize = config.systemSize;
    mechanism = config.mechanism;
    reltol = config.reltol;
    abstol = config.abstol;
    pressure = config.pressure;

}

void CanteraIntegrator::integrate(double t0, double t) {

}

void CanteraIntegrator::integrateSystem(int i) {

}