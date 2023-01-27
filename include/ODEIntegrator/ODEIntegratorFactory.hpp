#pragma once

#include <string>
#include <memory>
#include <stdexcept>

#include "yaml-cpp/yaml.h"

// Context and OutFileService
#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/Context/OutFileService.hpp"
// Loggers
#include "ODEIntegrator/Logger/Logger.hpp"
// Readers
#include "ODEIntegrator/InputOutput/Reader/Reader.hpp"
#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
// Integrators
#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Integrators/CanteraIntegrator.hpp"
#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"
#include "ODEIntegrator/Integrators/CVodeIntegratorOMP.hpp"
#include "ODEIntegrator/Integrators/CVodeIntegratorGPU.hpp"
// Writers
#include "ODEIntegrator/InputOutput/Writer/Writer.hpp"
#include "ODEIntegrator/InputOutput/Writer/csvWriter.hpp"

class ODEIntegratorFactory {
    public:
        ODEIntegratorFactory(YAML::Node _config): config(_config) {}

        std::shared_ptr<Context> createContext() const {
            std::shared_ptr<OutFileService> outFileService = createOutFileService();
            std::shared_ptr<BaseLogger> logger = createLogger(outFileService);
            return std::make_shared<Context>(outFileService, logger);
        }

        std::unique_ptr<Reader> createReader(std::shared_ptr<Context> ctx) const {
            std::string readerType = config["reader"]["type"].as<std::string>();
            if (readerType == "csvReader") {
                std::string inputCsvFile = config["reader"]["filename"].as<std::string>();
                return std::make_unique<csvReader>(ctx, inputCsvFile);
            }
            throw std::runtime_error("Reader type is not valid");
        }

        std::unique_ptr<Writer> createWriter(std::shared_ptr<Context> ctx) const {
            std::string writerType = config["writer"]["type"].as<std::string>();
            if (writerType == "csvWriter") {
                std::string outputCsvFile = config["writer"]["filename"].as<std::string>();
                return std::make_unique<csvWriter>(ctx, outputCsvFile);
            }
            throw std::runtime_error("Writer type is not valid");
        }

        std::unique_ptr<Integrator> createIntegrator() const {
            std::string integratorType = config["integrator"]["type"].as<std::string>();
            if (integratorType == "CanteraIntegrator") {
                return std::make_unique<CanteraIntegrator>();
            } else if (integratorType == "CanteraIntegratorOMP") {
                return std::make_unique<CanteraIntegratorOMP>();
            } else if (integratorType == "CVodeIntegrator") {
                return std::make_unique<CVodeIntegrator>();
            } else if (integratorType == "CVodeIntegratorOMP") {
                return std::make_unique<CVodeIntegratorOMP>();
            } else if (integratorType == "CVodeIntegratorGPU") {
                return std::make_unique<CVodeIntegratorGPU>();
            }
            throw std::runtime_error("Integrator type is not valid");
        }


    private:
        YAML::Node config;

        std::shared_ptr<OutFileService> createOutFileService() const {
            YAML::Node outFolder = config["outFileService"]["outFolder"];
            if (outFolder.IsNull())
                return std::make_shared<OutFileService>();
            else
                return std::make_shared<OutFileService>(outFolder.as<std::string>());
        }

        std::shared_ptr<BaseLogger> createLogger(std::shared_ptr<OutFileService> outFileService) const {
            // Unorder map to match string with enum class
            std::unordered_map<std::string, LogLevel> enumMap = {
                {"DEBUG", LogLevel::DEBUG},
                {"INFO", LogLevel::INFO}
            };
            // Set logLevel based on the config info
            LogLevel logLevel = enumMap[config["logger"]["logLevel"].as<std::string>()];

            std::string loggerType = config["logger"]["type"].as<std::string>();
            if (loggerType == "FileLogger")
                return std::make_shared<FileLogger>(logLevel, outFileService);
            else if (loggerType == "ConsoleLogger")
                return std::make_shared<ConsoleLogger>(logLevel);
            
            throw std::runtime_error("Logger type is not valid");
        }
};