#pragma once

#include <string>
#include <filesystem>
#include <fstream>

#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/Logger/BaseLogger.hpp"

#define info(x) info(x, __LINE__, __FILE__)
#define debug(x) debug(x, __LINE__, __FILE__)
#define error(x) error(x, __LINE__, __FILE__)

class FileLogger : public BaseLogger {
    protected:
        void log(std::string data) override;

    public:
        // TODO: remove context from logger
        FileLogger(LogLevel _logLevel, Context ctx);
        ~FileLogger();

    private:

        std::ofstream fout;
};

class ConsoleLogger : public BaseLogger {
    protected:
        void log(std::string data) override;

    public:
        ConsoleLogger(LogLevel _logLevel);
        ~ConsoleLogger();
};
