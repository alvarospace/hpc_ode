#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Integrators/CanteraIntegrator.hpp"


void CanteraIntegrator::init(IntegratorConfig config) {
    totalSize = config.totalSize;
    systemSize = config.systemSize;
    mechanism = config.mechanism;
    reltol = config.reltol;
    abstol = config.abstol;
    pressure = config.pressure;

    
}