#pragma once

#include "ODEIntegrator/Integrators/Integrator.hpp"

class CVodeIntegrator : public Integrator {
    public:
        virtual void dydt() = 0;
        virtual void jacobian() = 0;
};