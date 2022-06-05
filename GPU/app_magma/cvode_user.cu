#include <app_magma/cvode_user.cuh>

int eval_jacob_cvode(double t, N_Vector y, N_Vector ydot, SUNMatrix J, void* userData, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData *uData = (UserData*) userData;
  realtype *Jptr, *yptr, *JptrPy, *yptrPy;

  Jptr = SUNMatrix_MagmaDense_Data(J);
  yptr = N_VGetDeviceArrayPointer(y);

  /* Length of batched Jacobian matrix (should be = NSP * NSP * gpu_points) */
  // sunindextype JLength = SUNMatrix_MagmaDense_LData(J);

  JptrPy = uData->h_mem->jac;
  yptrPy = uData->h_mem->y;

  
  // Each GPU thread evaluate 1 jacobian matrix (1 system per thread)
  size_t nBlocks = (int) ceil( ((float) uData->nSystems) / BLOCKSIZE );
  dim3 dimGrid ( nBlocks );
  dim3 dimBlock ( BLOCKSIZE );

  /* Kernel call */
  kernel_eval_jacob<<< dimGrid, dimBlock >>>(uData->nSystems, t, uData->Pressure, yptr, Jptr, yptrPy, JptrPy, uData->d_mem);
  
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
  realtype *yptr, *ydotptr, *yptrPy, *ydotptrPy;

  yptr = N_VGetDeviceArrayPointer(y);
  ydotptr = N_VGetDeviceArrayPointer(ydot);

  yptrPy = uData->h_mem->y;
  ydotptrPy = uData->h_mem->dy;

  // Each GPU thread evaluate 1 dydt system
  size_t nBlocks = (int) ceil( ((float) uData->nSystems) / BLOCKSIZE );
  dim3 dimGrid ( nBlocks );
  dim3 dimBlock ( BLOCKSIZE );

  /* Kernel Call */
  kernel_dydt<<< dimGrid, dimBlock >>>(uData->nSystems, t, uData->Pressure, yptr, ydotptr, yptrPy, ydotptrPy, uData->d_mem);
  
  cudaDeviceSynchronize();
  cudaError_t cudaErr = cudaGetLastError();
  if (cudaErr != cudaSuccess) {
    fprintf(stderr, "\t ERROR in 'dydt_cvode': cudaGetLastError returned %s", cudaGetErrorName(cudaErr));
    return -1;
  }

  return 0;
}

__global__ void kernel_dydt(const int nSystems, const double t, const double P, double *ySun, double* dySun,
                            double *yPy, double *dyPy, const mechanism_memory *d_mem) {
  if (T_ID < nSystems) {
    // Reorder data for PyJac
    sun_to_pyjac_Y(ySun, yPy);

    dydt(t, P, yPy, dyPy, d_mem);

    // Reorder data back to Sundials
    pyjac_to_sun_Y(dyPy, dySun);
  }
}

__global__ void kernel_eval_jacob(const int nSystems, const double t, const double P, double *ySun, double *JSun,
                                  double *yPy, double *JPy, const mechanism_memory *d_mem) {

  if (T_ID < nSystems) {
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
  // yPy  = {T0,  T1,  T2, ...  TNSP, Y00, Y10, Y20, ..., Y(nSystems)0, Y01, Y11, Y21, ..., Y(nSystems)1, ...}

  // JSun = {J0, J1, ..., J(nSystems)} Each matrix is ordered column-major
  // JPy  = same structure that yPy, ordered column major

  for (int j = 0; j < NSP; j++) {
    // Location of the first system element for the current thread for "y" vector
    int sunSystemY = threadID * NSP;

    yPy[INDEX(j)] = ySun[sunSystemY + j];


    // Location of the first element of the system Jacobian matrix
    int sunSystemJac = threadID * NSP * NSP;

    for (int i = 0; i < NSP; i++) {
      JPy[INDEX(j*NSP + i)] =  JSun[sunSystemJac + j*NSP + i];
    }
  }
}

__device__ void pyjac_to_sun_YJ(double *yPy, double *ySun, double *JPy, double *JSun) {

  int threadID = threadIdx.x + blockIdx.x * blockDim.x;

    for (int j = 0; j < NSP; j++) {
    // Location of the first system element for the current thread for "y" vector
    int sunSystemY = threadID * NSP;

    ySun[sunSystemY + j] = yPy[INDEX(j)];

    // Location of the first element of the system Jacobian matrix
    int sunSystemJac = threadID * NSP * NSP;

    for (int i = 0; i < NSP; i++) {
      JSun[sunSystemJac + j*NSP + i] = JPy[INDEX(j * NSP + i)];
    }
  }

}

__device__ void sun_to_pyjac_Y(double *ySun, double *yPy) {

  int threadID = threadIdx.x + blockIdx.x * blockDim.x;

  // PyJac index -> #define INDEX(i) (T_ID + (i) * GRID_DIM)

  // ySun = {T0, Y00, Y01, ... Y0(NSP-1), T1, Y10, Y11, ... Y1(NSP-1), ...}
  // yPy  = {T0,  T1,  T2, ...  TNSP, Y00, Y10, Y20, ..., Y(nSystems)0, Y01, Y11, Y21, ..., Y(nSystems)1, ...}

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