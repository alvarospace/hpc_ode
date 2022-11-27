#include <string>
#include <filesystem>
#include <thread>
#include <sstream>
#include <chrono>

#include "ODEIntegrator/Logger/BaseLogger.hpp"

using std::string;

// Start queue_worker thread
BaseLogger::BaseLogger() {
    queue_worker = std::thread(&BaseLogger::processQueue, this);
}

// Wait for queue_worker to finish
BaseLogger::~BaseLogger() {
    log_queue.push(EXIT_FLAG);
    queue_worker.join();
}

void BaseLogger::buildMsg(string type, string data, int line, string file) {
    auto path = std::filesystem::path(file);
    auto thread_id = std::this_thread::get_id();
    std::stringstream msg_ss;
    msg_ss 
        << "[" << thread_id << "] " << getTimeStamp()
        << " [" << path.filename().string() << ":" << line << "] " 
        << type << data;
    log_queue.push(msg_ss.str());
}

string BaseLogger::getTimeStamp() const {
    // Time stamp
    const auto now = std::chrono::system_clock::now();
    const auto nowAsTimeT = std::chrono::system_clock::to_time_t(now);
    // Last miliseconds
    const auto nowMs = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    // Append the milliseconds to the time stamp(best precision in seconds)
    std::stringstream nowSs;
    nowSs
        << "[" << std::put_time(std::localtime(&nowAsTimeT), DFT_FMT.c_str())
        << ',' << std::setfill('0') << std::setw(3) << nowMs.count() << "]";
    return nowSs.str();
}

void BaseLogger::processQueue() {
    string msg;
    while (true) {
        log_queue.pop(msg);
        if (msg == EXIT_FLAG)
            break;
        log(msg);
    }
}
