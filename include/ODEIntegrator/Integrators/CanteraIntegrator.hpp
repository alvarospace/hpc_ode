#pragma once

#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

class CanteraIntegrator : public Integrator {
    public:
        void init(IntegratorConfig config) override;
        void integrate(double t0, double t) override;
        void clean() override {}

    private:
        void integrateSystem(double& temperature, double* species, double dt);
        void integrateSystem(double& temperature, double& enthalpy, double* species, double dt);
};