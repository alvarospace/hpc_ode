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
        virtual void init(std::shared_ptr<Context> _ctx, IntegratorConfig config);
        virtual void integrate(double t0, double t);
        virtual void clean();

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