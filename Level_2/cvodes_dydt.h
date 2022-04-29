#include "header.h"
#include "sundials/sundials_nvector.h"
#include "nvector/nvector_serial.h"
int dydt_cvodes(double t, N_Vector y, N_Vector ydot, void* f);
