#pragma once

#include <iostream>
#include <chrono>
#include <string>
#include <sstream>
#include <thread>
#include <mutex>
#include <filesystem>

#include "threadsafe_queue.hpp"

using namespace std;
using namespace std::chrono;

// TODO: Log thread id and __FILE__ and __LINE__

class Logger {
    public:
        void info(string data, string line, string file) {
            msg_pipe("[INFO]:  ", data, line, file);
        }

        void debug(string data, string line, string file) {
            msg_pipe("[DEBUG]: ", data, line, file);
        }

        void error(string data, string line, string file) {
            msg_pipe("[ERROR]: ", data, line, file);
        }

        Logger();
        ~Logger();

        // Delete copy constructors
        Logger(const Logger&) = delete;
        Logger& operator=(const Logger&) = delete;

    private:
        // Constants
        static inline string const DFT_FMT = "%T";
        static inline string const EXIT_FLAG = "SIGNAL EXIT";

        // Member properties
        ThreadSafeQueue<string> log_queue;
        thread queue_worker;
        
        // Private methods
        string time_stamp();
        void process_queue();
        void msg_pipe(string type, string data, string line, string file);

    protected:
        virtual void log(string data) {
            cout << data << endl;
        }
};

void Logger::msg_pipe(string type, string data, string line, string file) {
    auto thread_id = this_thread::get_id();
    stringstream msg_ss;
    msg_ss 
        << "[Thread - " << thread_id << "] " << time_stamp()
        << type << data;
    log_queue.push(msg_ss.str());
}

// Start queue_worker thread
Logger::Logger() {
    queue_worker = thread(&Logger::process_queue, this);
}

// Wait for queue_worker to finish
Logger::~Logger() {
    log_queue.push(EXIT_FLAG);
    queue_worker.join();
}

string Logger::time_stamp() {
    // Time stamp
    const auto now = system_clock::now();
    const auto nowAsTimeT = system_clock::to_time_t(now);
    // Last miliseconds
    const auto nowMs = duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    // Append the milliseconds to the time stamp(best precision in seconds)
    std::stringstream nowSs;
    nowSs
        << "[" << std::put_time(std::localtime(&nowAsTimeT), DFT_FMT.c_str())
        << '.' << std::setfill('0') << std::setw(3) << nowMs.count() << "] ";
    return nowSs.str();
}

void Logger::process_queue() {
    string msg;
    while (true) {
        log_queue.wait_pop(msg);
        if (msg == EXIT_FLAG)
            break;
        log(msg);
    }
}
