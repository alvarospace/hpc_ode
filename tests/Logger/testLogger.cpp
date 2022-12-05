#include <memory>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <omp.h>

#include "ODEIntegrator/Logger/Logger.hpp"
#include "ODEIntegrator/Context/OutFileService.hpp"

namespace fs = std::filesystem;

void testFileLoggerConcurrency() {
    std::string const OUTDIR = "./test_out/";
    auto fileService = std::make_shared<OutFileService>(OUTDIR);
    auto logger = std::make_unique<FileLogger>(LogLevel::INFO, fileService);

    // Print 20 logs with OpenMP threads
    omp_set_num_threads(10);
    #pragma omp parallel for
    for (int i = 0; i < 20; i++) {
        int id = omp_get_thread_num();
        logger->info("Thread ID: " + std::to_string(id) + " - Iter: " + std::to_string(i));
    }

    // Check that the number of logs are 20
    std::ifstream log_file;
    int lines = 0;
    fs::path logPath(fileService->getExecutionFolder());
    logPath /= "out.log";
    log_file.open(logPath.string());
    if (!log_file.is_open()) {
        throw std::runtime_error("Unable to open log file");
    }
    char buffer[100];
    while (log_file.getline(buffer, 100)) {
        lines++;
    }
    log_file.close();
    fs::remove_all(logPath.parent_path());

    // Num lines must be 20
    assert(lines == 20);
}


int main() {
    testFileLoggerConcurrency();
    return 0;
}