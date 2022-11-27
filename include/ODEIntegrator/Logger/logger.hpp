#pragma once

#include <iostream>
#include <chrono>
#include <string>
#include <ctime>
#include <thread>
#include <mutex>

#include "threadsafe_queue.hpp"

using namespace std;
using namespace std::chrono;

// TODO: Log thread id and __FILE__ and __LINE__

class Logger {
    public:
        void info(string data) {
            string msg = time_stamp() + "[INFO]:  " + data;
            log_queue.push(msg);
        }

        void debug(string data) {
            string msg = time_stamp() + "[DEBUG]:  " + data;
            log_queue.push(msg);
        }

        void error(string data) {
            string msg = time_stamp() + "[ERROR]:  " + data;
            log_queue.push(msg);
        }

        Logger();
        ~Logger();

        // Delete copy constructors
        Logger(const Logger&) = delete;
        Logger& operator=(const Logger&) = delete;

    private:
        static inline string const DFT_FMT = "%F %T";
        static inline int const TIME_BUFFER_MAX = 30;

        ThreadSafeQueue<string> log_queue;
        thread queue_worker;
        bool exit_flag;
        mutex exit_mut;


        string time_stamp();

        void process_queue();
    protected:
        virtual void log(string data) {
            cout << data << endl;
        }
};

// Start queue_worker thread
Logger::Logger() {
    exit_flag = false;
    queue_worker = thread(&Logger::process_queue, this);
}

// Wait for queue_worker to finish
Logger::~Logger() {
    exit_mut.lock();
    exit_flag = true;
    exit_mut.unlock();
    if (queue_worker.joinable())
        queue_worker.join();
}

string Logger::time_stamp() {
    // Time stamp
    auto date = system_clock::to_time_t(system_clock::now());

    char buffer [TIME_BUFFER_MAX];
    strftime(buffer, TIME_BUFFER_MAX, DFT_FMT.c_str(), gmtime(&date));

    string time_stamp(buffer);
    time_stamp = "[" + time_stamp + "] ";
    return time_stamp;
}

void Logger::process_queue() {
    string msg;
    unique_lock lk(exit_mut);
    while (!exit_flag) {
        lk.unlock();
        log_queue.wait_pop(msg);
        log(msg);
        lk.lock();
    }
}
