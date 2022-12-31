#pragma once

#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

class CanteraIntegrator : public Integrator {
    public:
        virtual void init(std::shared_ptr<Context> _ctx, IntegratorConfig config) override;
        virtual void integrate(double t0, double t) override;

    protected:
        void integrateSystem(double& temperature, double* species, double dt);
        void integrateSystem(double& temperature, double& enthalpy, double* species, double dt);
};

class CanteraIntegratorOMP : public CanteraIntegrator {
    public:
        void init(std::shared_ptr<Context> _ctx, IntegratorConfig config) override;
        void integrate(double t0, double t) override;
};