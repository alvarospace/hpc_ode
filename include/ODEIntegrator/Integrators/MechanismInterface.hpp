#pragma once

class MechanismInterface {
    public:
        virtual void dydt() = 0;
        virtual void jacobian() = 0;
};