#include <string>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <chrono>
#include <stdexcept>

#include "ODEIntegrator/Logger/BaseLogger.hpp"
#include "ODEIntegrator/Logger/Logger.hpp"
#include "ODEIntegrator/Context/Context.hpp"

namespace fs = std::filesystem;

/**************** Console Logger ******************/

ConsoleLogger::ConsoleLogger(LogLevel _logLevel) : BaseLogger(_logLevel) {}

ConsoleLogger::~ConsoleLogger() {
    finishWorker();
}

void ConsoleLogger::log(std::string data) {
    std::cout << data << std::endl;
}

/**************************************************/

/**************** File Logger *********************/
FileLogger::FileLogger(LogLevel _logLevel, Context ctx) : BaseLogger(_logLevel) {
    if (!ctx.isFolderReady())
        ctx.setUpFolder();
    
    auto logPath = fs::path(ctx.getOutFolder());

    // Formated TimeStamp
    auto const nowTimeT = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now()
    );
    std::stringstream logFileNameDated;
    logFileNameDated
        << "out" << "_"
        << std::put_time(std::localtime(&nowTimeT), "%F_%T")
        << ".log";

    // Log file dated path
    auto targetLogFile = logPath / logFileNameDated.str();

    fout.open(targetLogFile.string(), std::ios::out | std::ios::trunc);
    if (!fout.is_open()) {
        throw std::runtime_error("Error opening file for logging");
    }
}

FileLogger::~FileLogger() {
    finishWorker();
    fout.close();
}

void FileLogger::log(std::string data) {
    fout << data << std::endl;
}

/**************************************************/