#pragma once

#include <string>
#include <memory>

#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Logger/Logger.hpp"
#include "ODEIntegrator/Context/Context.hpp"

struct IntegratorConfig {
    std::string mechanism;
    double reltol;
    double abstol;
    double pressure;
};

//Abstract class that every integrator should inherit
class Integrator {
    public:
        virtual void init(Context ctx, IntegratorConfig config) {
            logger = ctx.getLogger();
            mesh = ctx.getMesh();
            totalSize = mesh->totalSize();
            systemSize = mesh->numSpecies();

            mechanism = config.mechanism;
            reltol = config.reltol;
            abstol = config.abstol;
            pressure = config.pressure;
        }
        virtual void integrate(double t0, double t) = 0;
        virtual void clean() = 0;

    protected:
        int totalSize;
        int systemSize;
        std::shared_ptr<Mesh> mesh;
        std::shared_ptr<BaseLogger> logger;

        std::string mechanism;
        double reltol;
        double abstol;
        double pressure;
};