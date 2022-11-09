#pragma once

#include <iostream>
#include <chrono>
#include <string>
#include <ctime>

using namespace std;
using namespace std::chrono;

class Logger
{
private:
    static inline string const DFT_FMT = "%F %T";
    static inline int const TIME_BUFFER_MAX = 30;

    string time_stamp();
protected:
    virtual void log(string data) {
        cout << data << endl;
    }

public:
    void info(string data) {
        log(time_stamp() + "[INFO]:  " + data);
    }

    void debug(string data) {
        log(time_stamp() + "[DEBUG]: " + data);
    }

    void error(string data) {
        log(time_stamp() + "[ERROR]: " + data);
    }
};

string Logger::time_stamp() {
    // Time stamp
    auto date = system_clock::to_time_t(system_clock::now());

    char buffer [TIME_BUFFER_MAX];
    strftime(buffer, TIME_BUFFER_MAX, DFT_FMT.c_str(), gmtime(&date));

    string time_stamp(buffer);
    time_stamp = "[" + time_stamp + "] ";
    return time_stamp;
}
