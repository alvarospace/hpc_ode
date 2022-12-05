#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Mesh/Point.hpp"
#include "ODEIntegrator/Mesh/PointBuilder.hpp"

using namespace std;

csvReader::csvReader(std::shared_ptr<Context> _ctx, std::string _csvFilename)
: Reader(_ctx) {
    string endWith {".csv"};
    if (_csvFilename.find(endWith) == string::npos) {
        stringstream ss;
        ss << "Input file should end with '.csv', actual: ";
        ss << _csvFilename;
        logger->error(ss.str());
        throw invalid_argument(ss.str());
    }
    csvFilename = _csvFilename;
    logger->info("csvReader created");
}

void csvReader::read() {
    stringstream ss;

    // Open csv file
    ifstream file(csvFilename);
    if (!file.is_open()) {
        ss << "Failing opening input file: ";
        ss << csvFilename;
        logger->error(ss.str());
        throw runtime_error(ss.str());
        ss.clear();
    }

    ss << "Reading file: \"" << csvFilename << "\"";
    logger->info(ss.str());
    ss.clear();

    // First line is the header
    // necessary to inspect the expected data
    string line;
    getline(file, line);
    HeaderInfo headerInfo = inspectHeader(line);
    ss << "NSP: " << headerInfo.nsp
       << "Temperature: " << headerInfo.hasTemperature << "\n"
       << "Coordenates: " << headerInfo.hasCoords << "\n"
       << "Enthalpy: "    << headerInfo.hasEnthalpy << "\n";
    logger->debug(ss.str());
    ss.clear();

    // Read all data
    while ( getline(file, line) ) {
        Point newPoint = readPoint(line, headerInfo);
        mesh->addPoint(newPoint);
    }
    file.close();
    logger->info("Read completed");
}

csvReader::HeaderInfo csvReader::inspectHeader(string header) {
    HeaderInfo headerInfo;

    // Format of parameters
    string nspFormat {"CO"};
    string enthalpyFormat {"ENTHA"};
    string temperatureFormat {"TEMPE"};
    string coordinatesFormat {"Points"};

    vector<string> columns;
    string column;
    stringstream stringIterator(header);
    while ( getline(stringIterator, column, ',') ) {
        columns.push_back(column);
    }

    int nspNum = count_if(begin(columns), end(columns), [&nspFormat](string column){
        return column.find(nspFormat) != string::npos;
    });
    headerInfo.nsp = nspNum;

    auto it = find_if(begin(columns), end(columns), [&enthalpyFormat](string column){
        return column.find(enthalpyFormat) != string::npos;
    });
    if (it != end(columns)){
        headerInfo.hasEnthalpy = true;
    }

    it = find_if(begin(columns), end(columns), [&temperatureFormat](string column){
        return column.find(temperatureFormat) != string::npos;
    });
    if (it != end(columns)){
        headerInfo.hasTemperature = true;
    }
    
    it = find_if(begin(columns), end(columns), [&coordinatesFormat](string column){
        return column.find(coordinatesFormat) != string::npos;
    });
    if (it != end(columns)){
        headerInfo.hasCoords = true;
    }

    return headerInfo;
}

Point csvReader::readPoint(string line, csvReader::HeaderInfo headerInfo) {
    stringstream lineItems(line);
    string item;

    // Species
    vector<double> species (headerInfo.nsp, 0.0f);
    for_each(begin(species), end(species), [&lineItems, &item](double &sp){
        getline(lineItems, item, ',');
        sp = stod(item);
    });

    // Enthalpy
    double enthalpy {};
    if (headerInfo.hasEnthalpy) {
        getline(lineItems, item, ',');
        enthalpy = stod(item);
    }

    // Temperature
    if (!headerInfo.hasTemperature){
        throw runtime_error("Input file does not have temperature!!");
    }
    getline(lineItems, item, ',');
    double temperature = stod(item);

    // Coordinates
    Coords coords;
    if (headerInfo.hasCoords) {
        getline(lineItems, item, ',');
        coords.x = stod(item);
        getline(lineItems, item, ',');
        coords.y = stod(item);
        getline(lineItems, item, ',');
        coords.z = stod(item);
    }

    // Create point
    PointBuilder pointBuilder(headerInfo.nsp);
    pointBuilder.addSpecies(species)
                .addTemperature(temperature);
    if (headerInfo.hasEnthalpy)
        pointBuilder.addEnthalpy(enthalpy);
    if (headerInfo.hasCoords)
        pointBuilder.addCoordinates(coords);

    return pointBuilder;
}