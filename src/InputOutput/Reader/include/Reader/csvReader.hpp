#pragma once

#include "Reader/Reader.hpp"
#include "Mesh/Mesh.hpp"
#include <string>

using std::string;

class csvReader : public Reader {
    private:
        string csvFilename;

        struct HeaderInfo {
            int nsp {};
            bool hasTemperature {false};
            bool hasCoords {false};
            bool hasEnthalpy {false};
        };

        HeaderInfo inspectHeader(string header);

        Point readPoint(string line, HeaderInfo headerInfo);
    
    public:
        csvReader(string _csvFilename) : csvFilename(_csvFilename) {}
        
        void read() override;
};