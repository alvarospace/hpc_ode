#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <stdexcept>
#include <sstream>

#include "ODEIntegrator/Logger/BaseLogger.hpp"
#include "ODEIntegrator/Logger/Logger.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Context/OutFileService.hpp"

class Context {
    public:
        // Default context
        Context(std::shared_ptr<OutFileService> _fileService) {
            fileService = _fileService;
            logger = std::make_shared<ConsoleLogger>(LogLevel::INFO);
            mesh = std::make_shared<Mesh>();
        }

        Context(std::shared_ptr<OutFileService> _fileService, std::shared_ptr<Mesh> _mesh) {
            fileService = _fileService;
            logger = std::make_shared<ConsoleLogger>(LogLevel::INFO);
            mesh = _mesh;
        }

        Context(std::shared_ptr<OutFileService> _fileService, std::shared_ptr<BaseLogger> _logger) {
            fileService = _fileService;
            logger = _logger;
            mesh = std::make_shared<Mesh>();
        }

        Context(std::shared_ptr<OutFileService> _fileService, std::shared_ptr<BaseLogger> _logger, std::shared_ptr<Mesh> _mesh) {
            fileService = _fileService;
            logger = _logger;
            mesh = _mesh;
        }

        void setLogger(std::shared_ptr<BaseLogger> _logger) {
            logger = _logger;
        }

        void setMesh(std::shared_ptr<Mesh> _mesh) {
            mesh = _mesh;
        }

        std::shared_ptr<BaseLogger> getLogger() const {
            return logger;
        }

        std::shared_ptr<Mesh> getMesh() const {
            return mesh;
        }

        std::shared_ptr<OutFileService> getOutFileService() const {
            return fileService;
        }

        std::shared_ptr<OutFileService> fileService;

    private:
        std::shared_ptr<BaseLogger> logger;
        std::shared_ptr<Mesh> mesh;
};