#pragma once

#include <memory>

#include "ODEIntegrator/Context/Context.hpp"

struct UserData {
    std::shared_ptr<Context> ctx;
    double pressure;
};