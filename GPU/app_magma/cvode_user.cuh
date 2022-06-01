// SUNDIALS
#include <cvode/cvode.h>
#include <sunmatrix/sunmatrix_magmadense.h>
#include <sunlinsol/sunlinsol_magmadense.h>
#include <nvector/nvector_cuda.h>
#include <sunmemory/sunmemory_cuda.h>

// PYJAC
extern "C"{
#include <libpyjacGPU/dydt.cuh>
#include <libpyjacGPU/jacob.cuh>
#include <libpyjacGPU/mechanism.cuh>
}

// CUDA block size
#define BLOCKSIZE NSP

/* Struct to CVode User Data */
struct UserData {
    sunindextype nBlocks;
    sunindextype nEquations;
    realtype Pressure;
    mechanism_memory *h_mem;
    mechanism_memory *d_mem;
};


int eval_jacob_cvode(realtype t, N_Vector y, N_Vector ydot, SUNMatrix Jac, void* f, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int dydt_cvode(realtype t, N_Vector y, N_Vector ydot, void* f);
__global__ kernel_dydt(const int nEquations, const double t, const double P, const double *y, const double* dy, const mechanism_memory *d_mem);

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