#pragma once

#include <memory>

#include "nvector/nvector_cuda.h"
#include "sunmatrix/sunmatrix_magmadense.h"

#include "ODEIntegrator/Context/Context.hpp"
#include "Mechanism/GPU/gpu_memory.cuh"

int dydt_cvode_GPU(double t, N_Vector y, N_Vector ydot, void* userdata);

int jacobian_cvode_GPU(double t, N_Vector y, N_Vector ydot, SUNMatrix J, void* userdata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int calc_gpu_points(std::shared_ptr<Context> ctx, int total_points, int &real_calculated_points);

void init_memory_gpu(std::shared_ptr<Context> ctx, int num_systems, mechanism_memory *pyjac_mem);

void free_memory_gpu(std::shared_ptr<Context> ctx, mechanism_memory *pyjac_mem);
