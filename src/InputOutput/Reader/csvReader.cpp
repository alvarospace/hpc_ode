#include "Reader/csvReader.hpp"
#include "Mesh/Mesh.hpp"
#include "Mesh/Point.hpp"
#include "Mesh/PointBuilder.hpp"
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

using namespace std;

void csvReader::read() {
    // Open csv file
    ifstream file(csvFilename);
    if (!file.is_open()) {
        stringstream ss;
        ss << "Failing opening file: ";
        ss << csvFilename << endl;
        throw runtime_error(ss.str());
    }

    // First line is the header
    // necessary to inspect the expected data
    string line;
    getline(file, line);
    HeaderInfo headerInfo = inspectHeader(line);

    // Read all data
    Mesh &mesh = Mesh::get();
    while ( getline(file, line) ) {
        Point newPoint = readPoint(line, headerInfo);
        mesh.addPoint(newPoint);
    }
    file.close();
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