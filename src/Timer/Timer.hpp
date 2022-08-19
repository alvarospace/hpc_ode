#pragma once

#include <chrono>

class Timer {
    public:
        void tic();
        void toc();
        double getTime() const;

    private:
        std::chrono::steady_clock::time_point t1;
        std::chrono::steady_clock::time_point t2;
        std::chrono::duration<double> measured_time;
};