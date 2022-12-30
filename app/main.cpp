#include <string>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <stdexcept>

#include "yaml-cpp/yaml.h"

#include "ODEIntegrator/ODEIntegratorHeader.hpp"

using namespace std;

// TODO: Finish the factory
class ODEIntegratorFactory {
    public:
        ODEIntegratorFactory(YAML::Node _config): config(_config) {}

        shared_ptr<Context> createContext() const {
            shared_ptr<OutFileService> outFileService = createOutFileService();
            shared_ptr<BaseLogger> logger = createLogger(outFileService);
            return make_shared<Context>(outFileService, logger);
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
    YAML::Node config = YAML::LoadFile("/home/almousa/TFM/hpc_cvode/app/config.yaml");

    ODEIntegratorFactory factory(config);

    shared_ptr<Context> ctx = factory.createContext();
    shared_ptr<BaseLogger> logger = ctx->getLogger();
    logger->info("Starting ODEApplication...");
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