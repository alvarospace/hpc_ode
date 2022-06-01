#include <app_magma/cvode_user.cuh>

int eval_jacob_cvode(double t, N_Vector y, N_Vector ydot, SUNMatrix jac, void* f, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{       
	double* local_y = NV_DATA_S(y);
  double* jacptr = SM_DATA_D(jac);
  
  /* Jacobian evaluation by PyJac is in column-major, same as Sundials */
	eval_jacob((double)t, *(double*)f, local_y, (double*)jacptr);

	return 0;
}

int dydt_cvode(realtype t, N_Vector y, N_Vector ydot, void* userData)
{
  UserData *uData;
  realtype *yptr, *ydotptr, *yptrPy, *ydotptrPy;

  uData = (UserData*) userData;
  yptrPy = uData->h_mem->y;
  ydotptrPy = uData->h_mem->dy;

  yptr = N_VGetDeviceArrayPointer(y);
  ydotptr = N_VGetDeviceArrayPointer(ydot);

  /* Copy initial condition to PyJac interface */
  cudaErrorCheck ( cudaMemcpy(yptrPy, yptr, uData->nEquations * sizeof(double), cudaMemcpyDeviceToDevice) );


  size_t nBlocks = (int) ceil( ((float) uData->nEquations) / BLOCKSIZE );
  dim3 dimGrid ( nBlocks );
  dim3 dimBlock ( BLOCKSIZE );

  /* Kernel Call */
  kernel_dydt<<< dimGrid, dimBlock >>>(uData->nEquations, t, uData->Pressure, yptrPy, ydotptrPy, uData->d_mem);
  
  cudaDeviceSynchronize();
  cudaError_t cudaErr = cudaGetLastError();
  if (cudaErr != cudaSuccess) {
    fprintf(stderr, "\t ERROR in 'dydt_cvode': cudaGetLastError returned %s", cudaGetErrorName(cudaErr));
    return -1;
  }

  /* Copy results to sundials interface */
  cudaErrorCheck ( cudaMemcpy(ydotptr, ydotptrPy, uData->nEquations * sizeof(double), cudaMemcpyDeviceToDevice) );

  return 0;
}

__global__ kernel_dydt(const int nEquations, const double t, const double P, const double *y, const double* dy, const mechanism_memory *d_mem) {
  if (T_ID < nEquations)
    dydt(t, P, y, dy, d_mem);
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