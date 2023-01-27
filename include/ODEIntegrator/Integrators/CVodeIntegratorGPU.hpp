#pragma once

#include <memory>

#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Integrators/MechanismDriver.hpp"
#include "ODEIntegrator/Integrators/Drivers/CVodeGPUDriver.cuh"
#include "ODEIntegrator/Integrators/CVodeGPUDataModels.hpp"

class CVodeIntegratorGPU : public Integrator, public MechanismDriver {
    public:
        void init(std::shared_ptr<Context> ctx, IntegratorConfig config) override;
        void integrate(double t0, double t) override;
        dydt_driver dydt_func() override;
        jacobian_driver jacobian_func() override;

    private:
        void integrateBatch(double dt, int processed_systems, int systems_to_process, int padded);
        void setGPUArrayInitialConditions(double *yptr, int processed_systems, int systems_to_process);
        void saveGPUArrayResults(double *yptr, int processed_systems, int systems_to_process);
        void check_return_value(void* returnValue, std::string const funcName, int const opt);

        std::unique_ptr<GPUUserData> uData;
};