#pragma once

#include <string>
#include <memory>

#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/InputOutput/Writer/Writer.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

class csvWriter : public Writer {
    public:
        csvWriter(std::shared_ptr<Context> _ctx, std::string _csvFilename);
        void write() override;

    private:
        std::string csvFilename;
};