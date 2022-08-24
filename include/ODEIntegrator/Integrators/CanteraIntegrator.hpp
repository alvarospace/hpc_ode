#pragma once

#include "ODEIntegrator/Integrators/Integrator.hpp"

class CanteraIntegrator : public Integrator {
    public:
        void init(IntegratorConfig config) override;
        void integrate(double t0, double t) override;
        void clean() override;
};