#pragma once

#include <string>
#include <memory>

#include "ODEIntegrator/InputOutput/Reader/Reader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Mesh/Point.hpp"


class csvReader : public Reader {
    private:
        std::string csvFilename;
        std::shared_ptr<Mesh> mesh;


        struct HeaderInfo {
            int nsp {};
            bool hasTemperature {false};
            bool hasCoords {false};
            bool hasEnthalpy {false};
        };

        HeaderInfo inspectHeader(std::string header);

        Point readPoint(std::string line, HeaderInfo headerInfo);
    
    public:
        csvReader(std::string _csvFilename, std::shared_ptr<Mesh> _mesh);
        
        void read() override;
};