#include <memory>
#include <string>
#include <fstream>

#include "ODEIntegrator/InputOutput/Writer/csvWriter.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Context/Context.hpp"

using namespace std;

csvWriter::csvWriter(shared_ptr<Context> _ctx, string _csvFilename)
: Writer(_ctx) {
    string endWith {".csv"};
    if (_csvFilename.find(endWith) == string::npos) {
        stringstream ss;
        ss << "Out file should end with '.csv', actual: ";
        ss << _csvFilename;
        throw invalid_argument(ss.str());
    }
    csvFilename = _csvFilename;
    logger->info("csvWriter created");
}