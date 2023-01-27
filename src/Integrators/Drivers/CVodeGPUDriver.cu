#include <memory>
#include <stdexcept>
#include <cmath>
#include <sstream>

#include "cuda_runtime.h"

#include "ODEIntegrator/Integrators/Drivers/CVodeGPUDriver.cuh"
#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/Integrators/CVodeGPUDataModels.hpp"

#include "Mechanism/GPU/gpu_memory.cuh"
#include "Mechanism/GPU/dydt.cuh"
#include "Mechanism/GPU/jacob.cuh"

// CUDA block size
#define BLOCKSIZE 32
// Max GPU memory allocation by PyJac
#define MAX_GPU_MEM_PYJAC 0.8

__global__ void kernel_dydt(const int nSystems, const double t, const double P, double *ySun, double* dySun, double *yPy, double *dyPy, mechanism_memory pyjac_mem);
__device__ void sun_to_pyjac_Y(double *ySun, double *yPy);
__device__ void pyjac_to_sun_Y(double *yPy, double *ySun);
__global__ void kernel_eval_jacob(const int nSystems, const double t, const double P, double *ySun, double *JSun, double *yPy, double *JPy, mechanism_memory pyjac_mem);
__device__ void sun_to_pyjac_YJ(double *ySun, double *yPy, double *JSun, double *JPy);
__device__ void pyjac_to_sun_YJ(double *yPy, double *ySun, double *JPy, double *JSun);

int dydt_cvode_GPU(realtype t, N_Vector y, N_Vector ydot, void* userData)
{
    GPUUserData* uData = static_cast<GPUUserData*>(userData);
    auto logger = uData->ctx->getLogger();
    realtype *yptr, *ydotptr, *yptrPy, *ydotptrPy;

    yptr = N_VGetDeviceArrayPointer(y);
    ydotptr = N_VGetDeviceArrayPointer(ydot);

    yptrPy = uData->pyjac_mem->y;
    ydotptrPy = uData->pyjac_mem->dy;

    mechanism_memory pyjac_mem = *(uData->pyjac_mem);

    // Each GPU thread process 1 dydt system
    size_t nBlocks = (int) ceil( ((float) uData->nSystems) / BLOCKSIZE );
    dim3 dimGrid ( nBlocks );
    dim3 dimBlock ( BLOCKSIZE );

    /* Kernel Call */
    kernel_dydt<<< dimGrid, dimBlock >>>(uData->nSystems, t, uData->pressure, yptr, ydotptr, yptrPy, ydotptrPy, pyjac_mem);
    cudaDeviceSynchronize();

    #ifdef TESTING
    //uData->test_y_sun_vs_py->ysun_vs_ypyjac();
    //uData->test_y_sun_vs_py->ysun_vs_dypyjac();
    #endif

    cudaError_t cudaErr = cudaGetLastError();
    if (cudaErr != cudaSuccess) {
        std::stringstream ss;
        ss << "cudaGetLastError returned: " << cudaGetErrorName(cudaErr);
        logger->error(ss.str());
        return -1;
    }

    return 0;
}

int jacobian_cvode_GPU(double t, N_Vector y, N_Vector ydot, SUNMatrix J, void* userData, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    GPUUserData* uData = static_cast<GPUUserData*>(userData);
    auto logger = uData->ctx->getLogger();
    realtype *Jptr, *yptr, *JptrPy, *yptrPy;

    Jptr = SUNMatrix_MagmaDense_Data(J);
    yptr = N_VGetDeviceArrayPointer_Cuda(y);

    /* Length of batched Jacobian matrix (should be = NSP * NSP * gpu_points) */
    // sunindextype JLength = SUNMatrix_MagmaDense_LData(J);

    JptrPy = uData->pyjac_mem->jac;
    yptrPy = uData->pyjac_mem->y;

    mechanism_memory pyjac_mem = *(uData->pyjac_mem);

    // Each GPU thread evaluate 1 jacobian matrix (1 system per thread)
    size_t nBlocks = (int) ceil( ((float) uData->nSystems) / BLOCKSIZE );
    dim3 dimGrid ( nBlocks );
    dim3 dimBlock ( BLOCKSIZE );

    /* Kernel call */
    kernel_eval_jacob<<< dimGrid, dimBlock >>>(uData->nSystems, t, uData->pressure, yptr, Jptr, yptrPy, JptrPy, pyjac_mem);

    #ifdef TESTING
    uData->test_jacobian->compare_matrices();
    #endif

    cudaDeviceSynchronize();
    cudaError_t cudaErr = cudaGetLastError();
    if (cudaErr != cudaSuccess) {
        std::stringstream ss;
        ss << "cudaGetLastError returned: " << cudaGetErrorName(cudaErr);
        logger->error(ss.str());
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
    std::stringstream ss;
    ss << "GPU systems allocated in this iteration: " << num_systems;
    logger->info(ss.str());
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

__global__ void kernel_dydt(const int nSystems, const double t, const double P, double *ySun, double* dySun,
                            double *yPy, double *dyPy, mechanism_memory pyjac_mem) {
    if (T_ID < nSystems) {
        // Give the device pointers to a device structure pointer
        mechanism_memory *d_mem = &pyjac_mem;

        // Reorder data for PyJac
        sun_to_pyjac_Y(ySun, yPy);

        dydt(t, P, yPy, dyPy, d_mem);

        // Reorder data back to Sundials
        pyjac_to_sun_Y(dyPy, dySun);
    }
}

__device__ void sun_to_pyjac_Y(double *ySun, double *yPy) {

    int threadID = threadIdx.x + blockIdx.x * blockDim.x;

    // PyJac index -> #define INDEX(i) (T_ID + (i) * GRID_DIM)

    // ySun = {T0, Y00, Y01, ... Y0(NSP-1), T1, Y10, Y11, ... Y1(NSP-1), ...}
    // yPy  = {T0,  T1,  T2, ...  T(nSystems), Y00, Y10, Y20, ..., Y(nSystems)0, Y01, Y11, Y21, ..., Y(nSystems)1, ...}

    for (int i = 0; i < NSP; i++) {
        yPy[INDEX(i)] = ySun[threadID * NSP + i];
    }
}

__device__ void pyjac_to_sun_Y(double *yPy, double *ySun) {

    int threadID = threadIdx.x + blockIdx.x * blockDim.x;

    for (int i = 0; i < NSP; i++) {
        ySun[threadID * NSP + i] = yPy[INDEX(i)];
    }
}

__global__ void kernel_eval_jacob(const int nSystems, const double t, const double P, double *ySun, double *JSun,
                                  double *yPy, double *JPy, mechanism_memory pyjac_mem) {

    if (T_ID < nSystems) {
        // Give the device pointers to a device structure pointer
        mechanism_memory *d_mem = &pyjac_mem;

        // Reorder data for PyJac
        sun_to_pyjac_YJ(ySun, yPy, JSun, JPy);

        // Jacobian analytic evaluation with PyJac
        eval_jacob(t, P, yPy, JPy, d_mem);

        // Reorder data back to Sundials
        pyjac_to_sun_YJ(yPy, ySun, JPy, JSun);
    }
}

__device__ void sun_to_pyjac_YJ(double *ySun, double *yPy, double *JSun, double *JPy) {

    int threadID = threadIdx.x + blockIdx.x * blockDim.x;

    // PyJac index -> #define INDEX(i) (T_ID + (i) * GRID_DIM)

    // ySun = {T0, Y00, Y01, ... Y0(NSP-1), T1, Y10, Y11, ... Y1(NSP-1), ...}
    // yPy  = {T0,  T1,  T2, ...  T(nSystems), Y00, Y10, Y20, ..., Y(nSystems)0, Y01, Y11, Y21, ..., Y(nSystems)1, ...}

    // JSun = {J0, J1, ..., J(nSystems)} Each matrix is ordered column-major
    // JPy  = same structure that yPy, ordered column major

    // Location of the first system element for the current thread for "y" vector
    int sunSystemY = threadID * NSP;

    // Location of the first element of the system Jacobian matrix
    int sunSystemJac = threadID * NSP * NSP;

    for (int j = 0; j < NSP; j++) {
        yPy[INDEX(j)] = ySun[sunSystemY + j];
        for (int i = 0; i < NSP; i++) {
            JPy[INDEX(j*NSP + i)] =  JSun[sunSystemJac + j*NSP + i];
        }
    }
}

__device__ void pyjac_to_sun_YJ(double *yPy, double *ySun, double *JPy, double *JSun) {

    int threadID = threadIdx.x + blockIdx.x * blockDim.x;

    // Location of the first system element for the current thread for "y" vector
    int sunSystemY = threadID * NSP;

    // Location of the first element of the system Jacobian matrix
    int sunSystemJac = threadID * NSP * NSP;

    for (int j = 0; j < NSP; j++) {
        ySun[sunSystemY + j] = yPy[INDEX(j)];
        for (int i = 0; i < NSP; i++) {
            JSun[sunSystemJac + j*NSP + i] = JPy[INDEX(j * NSP + i)];
        }
    }
}
