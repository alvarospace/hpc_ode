#include <memory>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <omp.h>

#include "ODEIntegrator/Logger/Logger.hpp"

namespace fs = std::filesystem;

void testFileLoggerConcurrency() {
    std::string const OUTDIR = "./test_out/";
    auto logger = std::make_unique<FileLogger>(LogLevel::INFO, OUTDIR);

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
    for (auto const& dir_entry : fs::directory_iterator(OUTDIR)) {
        log_file.open(dir_entry.path().string());
        if (!log_file.is_open()) {
            throw std::runtime_error("Unable to open log file");
        }
        char buffer[100];
        while (log_file.getline(buffer, 100)) {
            lines++;
        }
        log_file.close();
    }
    fs::remove_all(OUTDIR);

    // Num lines must be 20
    assert(lines == 20);
}


int main() {
    testFileLoggerConcurrency();
    return 0;
}