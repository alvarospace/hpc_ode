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

class Logger {
    public:
        void info(string data, int line, string file) {
            buildMsg("INFO:  ", data, line, file);
        }

        void debug(string data, int line, string file) {
            buildMsg("DEBUG: ", data, line, file);
        }

        void error(string data, int line, string file) {
            buildMsg("ERROR: ", data, line, file);
        }

        Logger();
        ~Logger();

        // Delete copy constructors
        Logger(const Logger&) = delete;
        Logger& operator=(const Logger&) = delete;

    private:
        // Constants
        static inline string const DFT_FMT = "%F %T";
        static inline string const EXIT_FLAG = "SIGNAL EXIT";

        // Member properties
        ThreadSafeQueue<string> log_queue;
        thread queue_worker;
        
        // Private methods
        string getTimeStamp();
        void processQueue();
        void buildMsg(string type, string data, int line, string file);

    protected:
        virtual void log(string data) {
            cout << data << endl;
        }
};

void Logger::buildMsg(string type, string data, int line, string file) {
    auto path = filesystem::path(file);
    auto thread_id = this_thread::get_id();
    stringstream msg_ss;
    msg_ss 
        << "[" << thread_id << "] " << getTimeStamp()
        << " [" << path.filename().string() << ":" << line << "] " 
        << type << data;
    log_queue.push(msg_ss.str());
}

// Start queue_worker thread
Logger::Logger() {
    queue_worker = thread(&Logger::processQueue, this);
}

// Wait for queue_worker to finish
Logger::~Logger() {
    log_queue.push(EXIT_FLAG);
    queue_worker.join();
}

string Logger::getTimeStamp() {
    // Time stamp
    const auto now = system_clock::now();
    const auto nowAsTimeT = system_clock::to_time_t(now);
    // Last miliseconds
    const auto nowMs = duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    // Append the milliseconds to the time stamp(best precision in seconds)
    std::stringstream nowSs;
    nowSs
        << "[" << std::put_time(std::localtime(&nowAsTimeT), DFT_FMT.c_str())
        << ',' << std::setfill('0') << std::setw(3) << nowMs.count() << "]";
    return nowSs.str();
}

void Logger::processQueue() {
    string msg;
    while (true) {
        log_queue.pop(msg);
        if (msg == EXIT_FLAG)
            break;
        log(msg);
    }
}
