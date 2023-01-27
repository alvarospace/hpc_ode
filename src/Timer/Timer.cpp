#include <chrono>

#include "ODEIntegrator/Timer/Timer.hpp"

void Timer::tic() {
    t1 = std::chrono::steady_clock::now();
}

void Timer::toc() {
    t2 = std::chrono::steady_clock::now();
    measured_time = t2 - t1;
}

double Timer::getTime() const {
    return measured_time.count();
}