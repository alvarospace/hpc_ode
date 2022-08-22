#include "Reader/csvReader.hpp"
#include "Mesh/Mesh.hpp"
#include <cassert>
#include <stdexcept>
#include <string>

using std::string;

void testReadGoodFile() {
    csvReader reader("data/good_input_file.csv");
    reader.read();
    Mesh &mesh = Mesh::get();
    assert(mesh.totalSize() > 0);
}

void testReadBadFile() {
    string const expectedMsg = "Input file does not have temperature!!";
    csvReader reader("data/file_without_temperature.csv");
    try {
        reader.read();
    } catch (std::runtime_error const &e) {
        assert(e.what() == expectedMsg);
    }
}

void testReadNotExistingFile() {
    string const expectedMsg = "Failing opening file: data/nothing.csv";
    csvReader reader("data/nothing.csv");
    try {
        reader.read();
    } catch (std::runtime_error const &e) {
        assert(e.what() == expectedMsg);
    }
}

void testReadTxtFile() {
    string const expectedMsg = "Input file should end with '.csv', actual: data/dummmy.txt";
    try {
        csvReader reader("data/dummmy.txt");
    } catch (std::invalid_argument const &e) {
        assert(e.what() == expectedMsg);
    }
}

int main() {
    testReadGoodFile();
    testReadBadFile();
    testReadNotExistingFile();
    testReadTxtFile();
    return 0;
}