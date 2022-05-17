// SUNDIALS
#include <cvode/cvode.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>

// PYJAC
extern "C"{
#include <libpyjac/dydt.h>
#include <libpyjac/jacob.h>
}

int eval_jacob_cvode(realtype t, N_Vector y, N_Vector ydot, SUNMatrix Jac, void* f, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int dydt_cvode(realtype t, N_Vector y, N_Vector ydot, void* f);

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