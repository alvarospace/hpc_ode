#pragma once

#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"

int dydt_cvode_serial(double t, N_Vector y, N_Vector ydot, void* userdata);

int jacobian_cvode_serial(double t, N_Vector y, N_Vector ydot, SUNMatrix J, void* userdata, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
