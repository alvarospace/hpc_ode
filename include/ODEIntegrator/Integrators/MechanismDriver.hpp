#pragma once

#include "sundials/sundials_nvector.h"
#include "sundials/sundials_matrix.h"

typedef int (*dydt_driver)(double t, N_Vector y, N_Vector ydot, void* user_data);

typedef int (*jacobian_driver)(double t, N_Vector y, N_Vector fy, SUNMatrix Jac,
             void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

class MechanismDriver {
    public:
        virtual dydt_driver dydt_func() = 0;
        virtual jacobian_driver jacobian_func() = 0;
};