#pragma once

#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/Logger/BaseLogger.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

// Abstract class that every writer should inherit
class Writer {
    public:
        // Pure virtual read function
        virtual void write() = 0;

        Writer(std::shared_ptr<Context> _ctx) {
            ctx = _ctx;
            mesh = ctx->getMesh();
            logger = ctx->getLogger();
        }

    protected:
        std::shared_ptr<Context> ctx;
        std::shared_ptr<Mesh> mesh;
        std::shared_ptr<BaseLogger> logger;
};