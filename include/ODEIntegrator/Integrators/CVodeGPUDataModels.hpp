#pragma once

#include <memory>

#include "ODEIntegrator/Context/Context.hpp"
#include "Mechanism/GPU/mechanism.cuh"
#include "sundials/sundials_types.h"

struct GPUUserData {
    std::shared_ptr<Context> ctx;
    double pressure;
    std::unique_ptr<mechanism_memory> pyjac_mem;
    sunindextype nSystems;
    sunindextype nEquations;
};