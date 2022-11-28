#include <string>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <chrono>
#include <stdexcept>

#include "ODEIntegrator/Logger/BaseLogger.hpp"
#include "ODEIntegrator/Logger/Logger.hpp"

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
FileLogger::FileLogger(LogLevel _logLevel, std::string logFileName) : BaseLogger(_logLevel) {
    auto logPath = fs::path(logFileName);
    auto parentPath = logPath.parent_path();

    // Create directory if it does not exist
    if (fs::exists(parentPath)){
        // Remove if the path exists and it's not a directory
        if (!fs::is_directory(parentPath)) {
            fs::remove_all(parentPath);
            fs::create_directory(parentPath);
        }
    } else {
        fs::create_directories(parentPath);
    }

    // Formated TimeStamp
    auto const nowTimeT = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now()
    );
    std::stringstream logFileNameDated;
    logFileNameDated
        << logPath.stem().string() << "_"
        << std::put_time(std::localtime(&nowTimeT), "%F_%T")
        << logPath.extension().string();

    // Log file dated path
    auto targetLogFile = fs::path(parentPath) / logFileNameDated.str();

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