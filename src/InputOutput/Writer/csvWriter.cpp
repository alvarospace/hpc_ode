#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "ODEIntegrator/InputOutput/Writer/csvWriter.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Context/Context.hpp"

using namespace std;

csvWriter::csvWriter(shared_ptr<Context> _ctx, string _csvFilename)
: Writer(_ctx) {
    string endWith {".csv"};
    if (_csvFilename.find(endWith) == string::npos) {
        stringstream ss;
        ss << "Output file should end with '.csv', actual: ";
        ss << _csvFilename;
        logger->error(ss.str());
        throw invalid_argument(ss.str());
    }
    csvFilename = _csvFilename;
    logger->info("csvWriter created");
}

void csvWriter::write() {
    stringstream ss;

    // Open file in write mode
    ofstream file(csvFilename);
    if (!file.is_open()) {
        ss << "Failing opening output file: ";
        ss << csvFilename;
        logger->error(ss.str());
        throw runtime_error(ss.str());
        ss.clear();
    }

    ss << "Writing file: \"" << csvFilename << "\"";
    logger->info(ss.str());
    ss.clear();

    // Write header
    writeHeader(file);

    // Write data
    int points  = mesh->totalSize();
    for (int i = 0; i < points; i++) {
        writePoint(file, i);
    }

    file.close();
    logger->info("Write completed");
}

void csvWriter::writePoint(ofstream& file, int const index) {
    file << scientific;
    // Species
    int nsp = mesh->numSpecies();
    vector species = mesh->getSpeciesVector(index);
    for (int i = 0; i < nsp; i++) {
        file << species[i] << ",";
    }

    // Enthalpy
    if (mesh->hasEnthalpy())
        file << mesh->getEnthalpyVector()[index] << ",";
    
    // Temperature
    file << defaultfloat << mesh->getTemperatureVector()[index];

    // Coordinates
    if (mesh->hasCoordinates()) {
        Coords coords = mesh->getCoordinatesVector()[index];
        file << "," << coords.x << "," << coords.y << "," << coords.z;
    }

    file << endl;
}

void csvWriter::writeHeader(ofstream& file) {
    file << setw(3) << setfill('0');

    // Species
    int nsp = mesh->numSpecies();
    for (int i = 0; i < nsp; i++) {
        file << "CO" << i << ",";
    }
    
    // Enthalpy
    if (mesh->hasEnthalpy())
        file << "ENTHA" << ",";
    
    // Temperature
    file << "TEMPE";

    // Coordinates
    if (mesh->hasCoordinates()) {
        file << setw(0);
        for (int i = 0; i < 3; i++) {
            file << "," << "Points:" << i;
        }
    }
    
    file << endl;
}