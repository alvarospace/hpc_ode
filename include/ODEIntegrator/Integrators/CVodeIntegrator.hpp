#pragma once

#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Integrators/MechanismDriver.hpp"

class CVodeIntegrator : public Integrator, public MechanismDriver {
    public:
        virtual void init(IntegratorConfig config) override;
        virtual void integrate(double t0, double t) override;
        virtual void clean() override;
        virtual void dydt_func() override;
        virtual void jacobian_func() override;
};