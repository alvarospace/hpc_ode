#include <cassert>
#include <stdexcept>
#include <string>
#include <memory>

#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

using std::string;

void testReadGoodFile() {
    auto mesh = std::make_shared<Mesh>();
    csvReader reader("data/good_input_file.csv", mesh);
    reader.read();
    assert(mesh->totalSize() > 0);
}

void testReadBadFile() {
    string const expectedMsg = "Input file does not have temperature!!";
    auto mesh = std::make_shared<Mesh>();
    csvReader reader("data/file_without_temperature.csv", mesh);
    try {
        reader.read();
    } catch (std::runtime_error const &e) {
        assert(e.what() == expectedMsg);
    }
}

void testReadNotExistingFile() {
    string const expectedMsg = "Failing opening file: data/nothing.csv";
    auto mesh = std::make_shared<Mesh>();
    csvReader reader("data/nothing.csv", mesh);
    try {
        reader.read();
    } catch (std::runtime_error const &e) {
        assert(e.what() == expectedMsg);
    }
}

void testReadTxtFile() {
    string const expectedMsg = "Input file should end with '.csv', actual: data/dummmy.txt";
    auto mesh = std::make_shared<Mesh>();
    try {
        csvReader reader("data/dummmy.txt", mesh);
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