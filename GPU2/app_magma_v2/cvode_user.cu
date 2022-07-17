#include <app_magma_v2/cvode_user.cuh>

int eval_jacob_cvode(double t, N_Vector y, N_Vector ydot, SUNMatrix J, void* userData, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData *uData = (UserData*) userData;
  
  // Each GPU thread evaluate 1 jacobian matrix (1 system per thread)
  size_t nBlocks = (int) ceil( ((float) uData->nSystems) / BLOCKSIZE );
  dim3 dimGrid ( nBlocks );
  dim3 dimBlock ( BLOCKSIZE );

  /* Kernel call */
  kernel_eval_jacob<<< dimGrid, dimBlock >>>(uData->nSystems, t, uData->Pressure, uData->d_mem->y, uData->d_mem->jac, uData->d_mem);
  
  cudaDeviceSynchronize();
  cudaError_t cudaErr = cudaGetLastError();
  if (cudaErr != cudaSuccess) {
    fprintf(stderr, "\t ERROR in 'dydt_cvode': cudaGetLastError returned %s", cudaGetErrorName(cudaErr));
    return -1;
  }

	return 0;
}

int dydt_cvode(realtype t, N_Vector y, N_Vector ydot, void* userData)
{
  UserData *uData = (UserData*) userData;

  // Each GPU thread evaluate 1 dydt system
  size_t nBlocks = (int) ceil( ((float) uData->nSystems) / BLOCKSIZE );
  dim3 dimGrid ( nBlocks );
  dim3 dimBlock ( BLOCKSIZE );

  /* Kernel Call */
  kernel_dydt<<< dimGrid, dimBlock >>>(uData->nSystems, t, uData->Pressure, uData->d_mem->y, uData->d_mem->dy, uData->d_mem);
  
  cudaDeviceSynchronize();
  cudaError_t cudaErr = cudaGetLastError();
  if (cudaErr != cudaSuccess) {
    fprintf(stderr, "\t ERROR in 'dydt_cvode': cudaGetLastError returned %s", cudaGetErrorName(cudaErr));
    return -1;
  }

  return 0;
}

__global__ void kernel_dydt(const int nSystems, const double t, const double P, double *y, double *dy, const mechanism_memory *d_mem) {
  if (T_ID < nSystems) 
    dydt(t, P, y, dy, d_mem);
}

__global__ void kernel_eval_jacob(const int nSystems, const double t, const double P, double *y, double *J, const mechanism_memory *d_mem) {

  if (T_ID < nSystems)
    eval_jacob(t, P, y, J, d_mem);
}




/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}