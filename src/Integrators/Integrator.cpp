#include <sstream>
#include <memory>

#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Logger/Logger.hpp"

void Integrator::init(std::shared_ptr<Context> _ctx, IntegratorConfig config) {
    ctx = _ctx;
    logger = ctx->getLogger();
    mesh = ctx->getMesh();
    totalSize = mesh->totalSize();
    systemSize = mesh->numSpecies();

    mechanism = config.mechanism;
    reltol = config.reltol;
    abstol = config.abstol;
    pressure = config.pressure;
    ompConfig = config.ompConfig;

    logger->info("Initializing integrator...");
    std::stringstream ss;
    ss << "{ mechanism: " << mechanism << ", "
        << "reltol: " << reltol << ", "
        << "abstol: " << abstol << ", "
        << "pressure: " << pressure << " }";
    logger->info(ss.str());
}

void Integrator::integrate(double t0, double t) {
    std::stringstream ss;
    ss << "Integrating interval of: " << t - t0 << " seconds";
    logger->info(ss.str());
}

void Integrator::clean() {
    logger->info("Integrator clean completed");
}