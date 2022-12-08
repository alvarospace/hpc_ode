// TODO: Start building GPU solver
#include <memory>

#include "ODEIntegrator/Integrators/MechanismDriver.hpp"
#include "ODEIntegrator/Integrators/CVodeIntegratorGPU.hpp"
#include "ODEIntegrator/Integrators/Drivers/CVodeGPUDriver.cuh"
#include "Mechanism/GPU/mechanism.cuh"

void CVodeIntegratorGPU::init(std::shared_ptr<Context> ctx, IntegratorConfig config) {
    // Do not call CVodeIntegrator init method but Integrator init
    Integrator::init(ctx, config);
    uData = std::make_unique<GPUUserData>();
    uData->pressure = pressure;
    uData->ctx = ctx;
    uData->pyjac_mem = std::make_unique<mechanism_memory>();
}


dydt_driver CVodeIntegratorGPU::dydt_func() {
    return dydt_cvode_GPU;
}

jacobian_driver CVodeIntegratorGPU::jacobian_func() {
    return jacobian_cvode_GPU;
}