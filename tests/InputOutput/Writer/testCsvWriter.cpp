#include <cassert>
#include <string>
#include <memory>
#include <fstream>
#include <filesystem>

#include "ODEIntegrator/Logger/Logger.hpp"
#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/InputOutput/Writer/csvWriter.hpp"

using namespace std;

void testCsvWriter() {
    auto fileService = make_shared<OutFileService>();
    auto logger = make_shared<FileLogger>(LogLevel::DEBUG, fileService);
    auto ctx = make_shared<Context>(fileService, logger);

    // Input
    csvReader reader(ctx, "./data/res_gri_10.csv");
    reader.read();

    // Output
    csvWriter writer(ctx, "gri_10_results.csv");
    writer.write();

    filesystem::path outPath(fileService->getExecutionFolder());
    outPath /= "gri_10_results.csv";
    ifstream file(outPath.string());
    string line;
    int lineCounter {0};
    while (getline(file, line)) {
        lineCounter++;
    }

    // Should be 10 systems + header
    assert(lineCounter == 11);
}

int main() {
    testCsvWriter();
    return 0;
}