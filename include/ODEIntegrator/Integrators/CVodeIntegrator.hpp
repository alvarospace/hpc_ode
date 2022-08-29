#pragma once

#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Integrators/MechanismInterface.hpp"

class CVodeIntegrator : public Integrator, public MechanismInterface {
    public:
        virtual void init(IntegratorConfig config) override;
        virtual void integrate(double t0, double t) override;
        virtual void clean() override;
        virtual void dydt() override;
        virtual void jacobian() override;
};