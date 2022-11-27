#pragma once

#include <string>
#include <thread>

#include "ODEIntegrator/Logger/ThreadSafeQueue.hpp"

class BaseLogger {
    protected:
        // Pure virtual function that implements the I/O
        virtual void log(std::string data) = 0;

    public:
        void info(std::string data, int line, std::string file) {
            buildMsg("INFO:  ", data, line, file);
        }

        void debug(std::string data, int line, std::string file) {
            buildMsg("DEBUG: ", data, line, file);
        }

        void error(std::string data, int line, std::string file) {
            buildMsg("ERROR: ", data, line, file);
        }

        BaseLogger();
        ~BaseLogger();

        // Delete copy constructors
        BaseLogger(const BaseLogger&) = delete;
        BaseLogger& operator=(const BaseLogger&) = delete;

    private:
        // Constants
        static inline std::string const DFT_FMT = "%F %T";
        static inline std::string const EXIT_FLAG = "SIGNAL EXIT";

        // Member properties
        ThreadSafeQueue<std::string> log_queue;
        std::thread queue_worker;
        
        // Private methods
        std::string getTimeStamp() const;
        void processQueue();
        void buildMsg(std::string type, std::string data, int line, std::string file);
};
