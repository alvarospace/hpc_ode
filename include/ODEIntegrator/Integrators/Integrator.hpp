#pragma once

#include <string>
#include <memory>
#include <sstream>

#include "yaml-cpp/yaml.h"

#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Logger/Logger.hpp"
#include "ODEIntegrator/Context/Context.hpp"

struct IntegratorConfig {
    std::string mechanism;
    double reltol;
    double abstol;
    double pressure;
    YAML::Node ompConfig;
};

//Abstract class that every integrator should inherit
class Integrator {
    public:
        virtual void init(std::shared_ptr<Context> _ctx, IntegratorConfig config) {
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
        virtual void integrate(double t0, double t) = 0;
        virtual void clean() = 0;

    protected:
        int totalSize;
        int systemSize;
        std::shared_ptr<Mesh> mesh;
        std::shared_ptr<BaseLogger> logger;
        std::shared_ptr<Context> ctx;

        std::string mechanism;
        double reltol;
        double abstol;
        double pressure;
        YAML::Node ompConfig;
};