#pragma once

#include <string>
#include <memory>

#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/InputOutput/Reader/Reader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Mesh/Point.hpp"

class csvReader : public Reader {
    private:
        std::string csvFilename;

        struct HeaderInfo {
            int nsp {};
            bool hasTemperature {false};
            bool hasCoords {false};
            bool hasEnthalpy {false};
        };

        HeaderInfo inspectHeader(std::string header);

        Point readPoint(std::string line, HeaderInfo headerInfo);
    
    public:
        csvReader(std::shared_ptr<Context> _ctx, std::string _csvFilename);
        
        void read() override;
};