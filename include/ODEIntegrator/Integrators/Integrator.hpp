#pragma once

#include <string>

struct IntegratorConfig {
    int totalSize;
    int systemSize;
    std::string mechanism;
    double reltol;
    double abstol;
    double pressure;
};


//Abstract class that every integrator should inherit
class Integrator {
    public:
        virtual void init(IntegratorConfig config) {
            totalSize = config.totalSize;
            systemSize = config.systemSize;
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
        std::string mechanism;
        double reltol;
        double abstol;
        double pressure;
};