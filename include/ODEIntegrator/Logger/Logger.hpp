#pragma once

#include <string>
#include <filesystem>
#include <fstream>

#include "ODEIntegrator/Logger/BaseLogger.hpp"
#include "ODEIntegrator/Context/OutFileService.hpp"

#define info(x) info(x, __LINE__, __FILE__)
#define debug(x) debug(x, __LINE__, __FILE__)
#define error(x) error(x, __LINE__, __FILE__)

class FileLogger : public BaseLogger {
    protected:
        void log(std::string data) override;

    public:
        FileLogger(LogLevel _logLevel, std::shared_ptr<OutFileService> fileService);
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
