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

    // // Print 20 logs with OpenMP threads
    // #pragma omp parallel for
    // for (int i = 0; i < 20; i++) {
    //     int id = omp_get_thread_num();
    //     logger->info("Thread ID: " + std::to_string(id) + " - Iter: " + std::to_string(i));
    // }

    // // Check that the number of logs are 20
    // for (auto const& dir_entry : fs::directory_iterator(OUTDIR)) {
    //     std::cout << dir_entry.path().string() << std::endl;
    // }

    // std::ifstream file("./test_out/");

}


int main() {
    testFileLoggerConcurrency();
    return 0;
}