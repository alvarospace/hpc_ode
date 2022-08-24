#pragma once

#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

class CanteraIntegrator : public Integrator {
    public:
        void init(IntegratorConfig config) override;
        void integrate(double t0, double t) override;
        void clean() override;

    private:
        void integrateSystem(int i);
        bool isCompatible();
};