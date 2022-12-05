#pragma once

#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"
#include "ODEIntegrator/Integrators/Drivers/CVodeGPUDriver.cuh"

class CVodeIntegratorGPU : public CVodeIntegrator {
    public:
        void init(std::shared_ptr<Context> _ctx, IntegratorConfig config) override;
        void integrate(double t0, double t) override;
        void clean() override;
        dydt_driver dydt_func() override;
        jacobian_driver jacobian_func() override;

    private:
        void integrateBatch();

        void initMemoryGPU(int num_systems);
};