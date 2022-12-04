#pragma once

#include <string>
#include <memory>

#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/InputOutput/Writer/Writer.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

class csvWriter : public Writer {
    public:
        csvWriter(std::string _csvFilename, Context _ctx);
        void write() override;

    private:
        Context ctx;
        std::string csvFilename;
        std::shared_ptr<Mesh> mesh;

        // std::string parsePoint()
};