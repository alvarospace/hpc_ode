#include <string>
#include <sstream>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <filesystem>
#include <fstream>

#include "yaml-cpp/yaml.h"

#include "ODEIntegrator/ODEIntegratorHeader.hpp"
#include "ODEIntegrator/ODEIntegratorFactory.hpp"

using namespace std;

void runODEApplication(string configFileName) {
    YAML::Node config = YAML::LoadFile(configFileName);
    
    // Timers
    Timer readTimer, writeTimer, integrateTimer;

    ODEIntegratorFactory factory(config);

    shared_ptr<Context> ctx = factory.createContext();
    shared_ptr<BaseLogger> logger = ctx->getLogger();
    logger->info("Starting ODEApplication...");

    // Save execution config.yaml
    shared_ptr<OutFileService> outFileService = ctx->getOutFileService();
    string outPath = outFileService->getExecutionFolder();
    filesystem::path fromFile(configFileName);
    filesystem::path toFile(outPath);
    toFile /= "config.yaml";
    logger->info("Copying config.yaml to results directory");
    filesystem::copy_file(fromFile, toFile);

    // Read data
    readTimer.tic();
    unique_ptr<Reader> reader = factory.createReader(ctx);
    reader->read();
    readTimer.toc();

    // Integrator
    integrateTimer.tic();
    IntegratorConfig integratorConfig;
    integratorConfig.reltol = config["integrator"]["reltol"].as<double>();
    integratorConfig.abstol = config["integrator"]["abstol"].as<double>();
    integratorConfig.mechanism = config["integrator"]["mechanism"].as<string>();
    integratorConfig.pressure = config["integrator"]["pressure"].as<double>();
    integratorConfig.ompConfig = config["integrator"]["omp"];
    double dt = config["integrator"]["dt"].as<double>();
    unique_ptr<Integrator> integrator = factory.createIntegrator();

    integrator->init(ctx, integratorConfig);
    integrator->integrate(0.0, dt);
    integrator->clean();
    integrateTimer.toc();

    // Writer
    writeTimer.tic();
    unique_ptr<Writer> writer = factory.createWriter(ctx);
    writer->write();
    writeTimer.toc();

    // Log and Write performance results
    stringstream ss;
    ss << "{ readTime: " << readTimer.getTime() << ", "
       << "integrationTime: " << integrateTimer.getTime() << ", "
       << "writeTime: " << writeTimer.getTime() << ", "
       << "totalTime: " 
       << readTimer.getTime() + integrateTimer.getTime() + writeTimer.getTime()
       << " }";
    logger->info(ss.str());
    logger->info("Writing performance.yaml to results directory");
    YAML::Node performanceYaml;
    performanceYaml["readTime"] = readTimer.getTime();
    performanceYaml["integrationTime"] = integrateTimer.getTime();
    performanceYaml["writeTime"] = writeTimer.getTime();
    performanceYaml["totalTime"] = readTimer.getTime() + integrateTimer.getTime() + writeTimer.getTime();
    
    filesystem::path performancePath(outPath);
    performancePath /= "performance.yaml";
    ofstream performanceFile(performancePath.string());
    performanceFile << performanceYaml;

    logger->info("End of ODEApplication execution");
}

void usage(string programName) {
    stringstream ss;
    ss
        << programName << " "
        << "<path to config.yaml>";
    cout << ss.str() << endl;
}

int main(int argc, char *argv[]) {
    string configFileName;
    // Command-line arguments
    if (argc == 2) {
        string arg1 = argv[1];
        if (arg1 == "-h" || arg1 == "--help") {
            usage(argv[0]);
            return 0;

        } else {
            configFileName = arg1;
        }
    } else {
        cout << "Error: No argument provided" << endl;
        return 1;
    }

    try {
        runODEApplication(configFileName);

    } catch (exception const &e) {
        e.what();
        return 1;
    }

    return 0;
}