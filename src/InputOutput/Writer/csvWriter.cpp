#include <memory>
#include <string>
#include <fstream>

#include "ODEIntegrator/InputOutput/Writer/csvWriter.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Context/Context.hpp"

using namespace std;

// TODO: Remove context from logger
csvWriter::csvWriter(string _csvFilename, Context _ctx) {
    ctx = _ctx;
    mesh = ctx.getMesh();

    string endWith {".csv"};
    if (_csvFilename.find(endWith) == string::npos) {
        stringstream ss;
        ss << "Out file should end with '.csv', actual: ";
        ss << _csvFilename;
        throw invalid_argument(ss.str());
    }
    csvFilename = _csvFilename;
}