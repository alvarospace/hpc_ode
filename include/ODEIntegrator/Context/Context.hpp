#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <stdexcept>
#include <sstream>

#include "ODEIntegrator/Logger/BaseLogger.hpp"
#include "ODEIntegrator/Logger/Logger.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

namespace fs = std::filesystem;

// TODO: Setup folder on constructor
class Context {
    public:
        // Default context
        Context() {
            logger = std::make_shared<ConsoleLogger>(LogLevel::INFO);
            mesh = std::make_shared<Mesh>();
        }

        Context(std::shared_ptr<Mesh> _mesh) {
            logger = std::make_shared<ConsoleLogger>(LogLevel::INFO);
            mesh = _mesh;
        }

        Context(std::shared_ptr<BaseLogger> _logger, std::shared_ptr<Mesh> _mesh) {
            logger = _logger;
            mesh = _mesh;
        }

        void setLogger(std::shared_ptr<BaseLogger> _logger) {
            logger = _logger;
        }

        void setMesh(std::shared_ptr<Mesh> _mesh) {
            mesh = _mesh;
        }

        std::shared_ptr<BaseLogger> getLogger() {
            return logger;
        }

        std::shared_ptr<Mesh> getMesh() {
            return mesh;
        }

        void setUpFolder() {
            // Create directory if it does not exist
            if (fs::exists(outFolder)){
                // Remove if the path exists and it's not a directory
                if (!fs::is_directory(outFolder)) {
                    fs::remove_all(outFolder);
                    fs::create_directory(outFolder);
                }
            } else {
                fs::create_directories(outFolder);
            }
            folderReady = true;
        }

        void setOutFolder(std::string _outFolder) {
            if (folderReady) 
                throw std::logic_error("Out folder already setup, cannot be changed");
            
            outFolder = _outFolder;
        }

        std::string getOutFolder() {
            return outFolder;
        }

        bool isFolderReady() {
            return folderReady;
        }
    
    private:
        std::shared_ptr<BaseLogger> logger;
        std::shared_ptr<Mesh> mesh;

        bool folderReady = false;

        std::string outFolder = "./out/";
};