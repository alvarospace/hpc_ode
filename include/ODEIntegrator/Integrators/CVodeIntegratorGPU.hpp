#pragma once

#include <memory>

#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"
#include "ODEIntegrator/Integrators/Drivers/CVodeGPUDriver.cuh"
#include "Mechanism/GPU/mechanism.cuh"

struct GPUUserData {
    std::shared_ptr<Context> ctx;
    double pressure;
    std::unique_ptr<mechanism_memory> pyjac_mem;
};

class CVodeIntegratorGPU : public CVodeIntegrator {
    public:
        void init(std::shared_ptr<Context> ctx, IntegratorConfig config) override;
        void integrate(double t0, double t) override;
        void clean() override;
        dydt_driver dydt_func() override;
        jacobian_driver jacobian_func() override;



    private:
        void integrateBatch();

        void initMemoryGPU(int num_systems);

        std::unique_ptr<GPUUserData> uData;
};