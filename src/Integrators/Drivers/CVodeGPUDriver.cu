#include <memory>
#include <stdexcept>
#include <cmath>

#include "cuda_runtime.h"

#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/Integrators/Drivers/CVodeGPUDriver.cuh"
#include "ODEIntegrator/Integrators/CVodeDataModels.hpp"
#include "Mechanism/GPU/gpu_memory.cuh"


// CUDA block size
#define BLOCKSIZE 32
// Max GPU memory allocation by PyJac
#define MAX_GPU_MEM_PYJAC 0.8

// TODO: Complete the GPU integrator and test
int dydt_cvode(realtype t, N_Vector y, N_Vector ydot, void* userData)
{
  GPUUserData* uData = (GPUUserData*) userData;
  realtype *yptr, *ydotptr, *yptrPy, *ydotptrPy;

  yptr = N_VGetDeviceArrayPointer(y);
  ydotptr = N_VGetDeviceArrayPointer(ydot);

  yptrPy = uData->pyjac_mem->y;
  ydotptrPy = uData->pyjac_mem->dy;

  mechanism_memory pyjac_mem = *(uData->pyjac_mem);

  // Each GPU thread evaluate 1 dydt system
  size_t nBlocks = (int) ceil( ((float) uData->nSystems) / BLOCKSIZE );
  dim3 dimGrid ( nBlocks );
  dim3 dimBlock ( BLOCKSIZE );

  /* Kernel Call */
  kernel_dydt<<< dimGrid, dimBlock >>>(uData->nSystems, t, uData->Pressure, yptr, ydotptr, yptrPy, ydotptrPy, pyjac_mem);
  cudaDeviceSynchronize();

  #ifdef TESTING
  //uData->test_y_sun_vs_py->ysun_vs_ypyjac();
  //uData->test_y_sun_vs_py->ysun_vs_dypyjac();
  #endif

  cudaError_t cudaErr = cudaGetLastError();
  if (cudaErr != cudaSuccess) {
    fprintf(stderr, "\t ERROR in 'dydt_cvode': cudaGetLastError returned %s", cudaGetErrorName(cudaErr));
    return -1;
  }

  return 0;
}


int calc_gpu_points(std::shared_ptr<Context> ctx, int total_points, int &real_calculated_points) {
    auto logger = ctx->getLogger();
    size_t mech_size = required_mechanism_size();
    size_t free_mem = 0;
    size_t total_mem = 0;

    cudaErrorCheck( cudaMemGetInfo(&free_mem, &total_mem) );

    int max_allocated_points = int(floor( MAX_GPU_MEM_PYJAC * ((double)free_mem / (double)mech_size) ));

    // Choose between the remaining points and the maximum allocatable 
    real_calculated_points = std::min(total_points, max_allocated_points);

    // Transform padded in a number multiple of BLOCKSIZE, ej: 1000 -> 1024
    int padded = int(std::ceil(real_calculated_points / float(BLOCKSIZE)) * BLOCKSIZE);

    if (padded == 0) {
        logger->error("Mechanism is too large, cannot allocate any point on GPU");
        throw std::runtime_error("Mechanism is too large, cannot allocate any point on GPU");
    }

    return padded;
}

void init_memory_gpu(std::shared_ptr<Context> ctx, int num_systems, mechanism_memory *pyjac_mem) {
    auto logger = ctx->getLogger();
    logger->info("Initializing PyJac GPU memory...");

    cudaError_t cuda_err;

    /* Allocate Memory */

    cuda_err = cudaMalloc((void **) &(pyjac_mem->y), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->dy), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->conc), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->fwd_rates), FWD_RATES * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->rev_rates), REV_RATES * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->spec_rates), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->cp), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->h), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->dBdT), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->jac), NSP * NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->var), 1 * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->J_nplusjplus), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->pres_mod), PRES_MOD_RATES * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMalloc error code: " + cuda_err);
    }

    /* Init memory */

    cuda_err = cudaMemset(pyjac_mem->spec_rates, 0, NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMemset error code: " + cuda_err);
    }

    cuda_err = cudaMemset(pyjac_mem->dy, 0, NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMemset error code: " + cuda_err);
    }

    cuda_err = cudaMemset(pyjac_mem->jac, 0, NSP * NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        logger->error("cudaMemset error code: " + cuda_err);
    }
    logger->info("GPU systems allocated in this iteration: " + num_systems);
    logger->info("PyJac GPU memory initialized");
}

void free_memory_gpu(std::shared_ptr<Context> ctx, mechanism_memory *pyjac_mem) {
    auto logger = ctx->getLogger();
    cudaError_t cuda_err;

    cuda_err = cudaFree(pyjac_mem->y);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->dy);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->conc);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->fwd_rates);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->rev_rates);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->spec_rates);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->cp);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->h);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->dBdT);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->jac);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->var);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->J_nplusjplus);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    cuda_err = cudaFree(pyjac_mem->pres_mod);
    if (cuda_err != cudaSuccess) {
        logger->error("cudaFree error code: " + cuda_err);
    }

    logger->info("PyJac GPU memory freed");
}
