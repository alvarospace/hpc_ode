// SUNDIALS
#include <cvode/cvode.h>
#include <sunmatrix/sunmatrix_magmadense.h>
#include <sunlinsol/sunlinsol_magmadense.h>
#include <nvector/nvector_cuda.h>
#include <sunmemory/sunmemory_cuda.h>
#include <cuda.h>
#include <cuda_runtime.h>

// PYJAC
#include <libpyjacGPU/dydt.cuh>
#include <libpyjacGPU/jacob.cuh>
#include <libpyjacGPU/mechanism.cuh>
#include <libpyjacGPU/gpu_macros.cuh>


// CUDA block size
#define BLOCKSIZE 16

/* Struct to CVode User Data */
struct UserData {
    sunindextype nSystems;
    sunindextype nEquations;
    realtype Pressure;
    mechanism_memory *h_mem;
    mechanism_memory *d_mem;
};


int eval_jacob_cvode(double t, N_Vector y, N_Vector ydot, SUNMatrix J, void* userData, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int dydt_cvode(realtype t, N_Vector y, N_Vector ydot, void* userData);
__global__ void kernel_dydt(const int nSystems, const double t, const double P, double *ySun, double* dySun, double *yPy, double *dyPy, const mechanism_memory *d_mem);
__global__ void kernel_eval_jacob(const int nSystems, const double t, const double P, double *ySun, double *JSun, double *yPy, double *JPy, const mechanism_memory *d_mem);
__device__ void sun_to_pyjac_YJ(double *ySun, double *yPy, double *JSun, double *JPy);
__device__ void pyjac_to_sun_YJ(double *yPy, double *ySun, double *JPy, double *JSun);
__device__ void sun_to_pyjac_Y(double *ySun, double *yPy);
__device__ void pyjac_to_sun_Y(double *yPy, double *ySun);


/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */
int check_retval(void *returnvalue, const char *funcname, int opt);