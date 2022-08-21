#include "Timer.hpp"
#include <cassert>
#include <unistd.h>

void testTimer() {
    int const sleepTime {1'000'000}; //Micro-seconds
    int const expectedTime = sleepTime / 1'000'000;

    Timer timer;
    timer.tic();
    usleep(sleepTime);
    timer.toc();

    int passedTime = static_cast<int>(timer.getTime());
    assert(passedTime == expectedTime);
}

int main() {
    testTimer();
    return 0;
}