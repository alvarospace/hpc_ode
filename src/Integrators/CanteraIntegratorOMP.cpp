#include "ODEIntegrator/Integrators/CanteraIntegrator.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

#include "Integrators/OpenMPRuntime.hpp"

void CanteraIntegratorOMP::init(std::shared_ptr<Context> _ctx, IntegratorConfig config) {
    CanteraIntegrator::init(_ctx, config);
    setOMPRuntime(config.ompConfig, logger);
    logger->info("CanteraIntegratorOMP initialized");
}

void CanteraIntegratorOMP::integrate(double t0, double t) {
    Integrator::integrate(t0, t);
    double* temperatures = nullptr;
    double* enthalpies = nullptr;
    if (mesh->hasEnthalpy())
        enthalpies = mesh->getEnthalpyPointer();

    temperatures = mesh->getTemperaturePointer();

    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < totalSize; i++) {
        double* species = mesh->getSpeciesPointer(i);
        if (mesh->hasEnthalpy())
            integrateSystem(temperatures[i], enthalpies[i], species, t);
        else
            integrateSystem(temperatures[i], species, t);
    }
}