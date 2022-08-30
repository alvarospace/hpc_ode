#pragma once

//typedef int (*dydt)(double t, N_Vector y, N_vector ydot, void* user_data);

class MechanismDriver {
    public:
        virtual void dydt_func() = 0;
        virtual void jacobian_func() = 0;
};