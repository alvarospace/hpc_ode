#include "csvReader.hpp"
#include "Mesh.hpp"
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

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
}