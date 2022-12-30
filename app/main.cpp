#include <string>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <filesystem>

#include "yaml-cpp/yaml.h"

#include "ODEIntegrator/ODEIntegratorHeader.hpp"

using namespace std;

class ODEIntegratorFactory {
    public:
        ODEIntegratorFactory(YAML::Node _config): config(_config) {}

        shared_ptr<Context> createContext() const {
            shared_ptr<OutFileService> outFileService = createOutFileService();
            shared_ptr<BaseLogger> logger = createLogger(outFileService);
            return make_shared<Context>(outFileService, logger);
        }

        unique_ptr<Reader> createReader(shared_ptr<Context> ctx) const {
            string readerType = config["reader"]["type"].as<string>();
            if (readerType == "csvReader") {
                string inputCsvFile = config["reader"]["filename"].as<string>();
                return make_unique<csvReader>(ctx, inputCsvFile);
            }
            throw runtime_error("Reader type is not valid");
        }

        unique_ptr<Writer> createWriter(shared_ptr<Context> ctx) const {
            string writerType = config["writer"]["type"].as<string>();
            if (writerType == "csvWriter") {
                string outputCsvFile = config["writer"]["filename"].as<string>();
                return make_unique<csvWriter>(ctx, outputCsvFile);
            }
            throw runtime_error("Writer type is not valid");
        }

        unique_ptr<Integrator> createIntegrator() const {
            string integratorType = config["integrator"]["type"].as<string>();
            if (integratorType == "CanteraIntegrator") {
                return make_unique<CanteraIntegrator>();
            } else if (integratorType == "CanteraIntegratorOMP") {
                return make_unique<CanteraIntegratorOMP>();
            } else if (integratorType == "CVodeIntegrator") {
                return make_unique<CVodeIntegrator>();
            } else if (integratorType == "CVodeIntegratorOMP") {
                return make_unique<CVodeIntegratorOMP>();
            } else if (integratorType == "CVodeIntegratorGPU") {
                return make_unique<CVodeIntegratorGPU>();
            }
            throw runtime_error("Integrator type is not valid");
        }


    private:
        YAML::Node config;

        shared_ptr<OutFileService> createOutFileService() const {
            YAML::Node outFolder = config["outFileService"]["outFolder"];
            if (outFolder.IsNull())
                return make_shared<OutFileService>();
            else
                return make_shared<OutFileService>(outFolder.as<string>());
        }

        shared_ptr<BaseLogger> createLogger(shared_ptr<OutFileService> outFileService) const {
            // Unorder map to match string with enum class
            unordered_map<string, LogLevel> enumMap = {
                {"DEBUG", LogLevel::DEBUG},
                {"INFO", LogLevel::INFO}
            };
            // Set logLevel based on the config info
            LogLevel logLevel = enumMap[config["logger"]["logLevel"].as<string>()];

            string loggerType = config["logger"]["type"].as<string>();
            if (loggerType == "FileLogger")
                return make_shared<FileLogger>(logLevel, outFileService);
            else if (loggerType == "ConsoleLogger")
                return make_shared<ConsoleLogger>(logLevel);
            
            throw runtime_error("Logger type is not valid");
        }
};

void runODEApplication(string configFileName) {
    YAML::Node config = YAML::LoadFile(configFileName);
    
    // Timers
    Timer readTimer, writeTimer, integrateTimer, totalTimer;

    // Start total timer
    totalTimer.tic();

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
    logger->info("Copying config.yaml to evidences");
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
    integrateTimer.toc();

    // Writer
    writeTimer.tic();
    unique_ptr<Writer> writer = factory.createWriter(ctx);
    writer->write();
    writeTimer.toc();

    integrator->clean();

    // End total timer
    totalTimer.toc();

    // TODO: Write the performance.yaml file with the times
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