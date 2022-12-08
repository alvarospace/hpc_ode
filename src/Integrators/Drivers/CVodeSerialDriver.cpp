#include "ODEIntegrator/Integrators/Drivers/CVodeSerialDriver.hpp"
#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"

#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"

// C code from mechanism
#ifdef __cplusplus
extern "C" {
#endif
#include "Mechanism/CPU/dydt.h"
#include "Mechanism/CPU/jacob.h"
#ifdef __cplusplus
}
#endif


int dydt_cvode_serial(double t, N_Vector y, N_Vector ydot, void* userdata) {
    double* yptr = N_VGetArrayPointer(y);
    double* ydotptr = N_VGetArrayPointer(ydot);
    UserData* udata = static_cast<UserData*>(userdata);

    dydt(t, udata->pressure, yptr, ydotptr);

    return 0;
}

int jacobian_cvode_serial(double t, N_Vector y, N_Vector ydot, SUNMatrix J, void* userdata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    double* yptr = N_VGetArrayPointer(y);
    double* Jptr = SM_DATA_D(J);
    UserData* udata = static_cast<UserData*>(userdata);

    eval_jacob(t, udata->pressure, yptr, Jptr);

    return 0;
}