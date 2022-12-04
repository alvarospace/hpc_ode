#pragma once

#include <memory>

#include "ODEIntegrator/Logger/Logger.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

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
    
    private:
        std::shared_ptr<BaseLogger> logger;
        std::shared_ptr<Mesh> mesh;
};