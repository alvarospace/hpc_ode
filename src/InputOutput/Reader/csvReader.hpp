#pragma once

#include "Reader.hpp"
#include "Mesh/Mesh.hpp"
#include <string>

using std::string;

class csvReader : public Reader {
    public:
        csvReader(string _csvFilename);
        void read() override;

    private:
        string csvFilename;
};