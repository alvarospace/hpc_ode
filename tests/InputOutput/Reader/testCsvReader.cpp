#include <cassert>
#include <stdexcept>
#include <string>
#include <memory>

#include "ODEIntegrator/Context/Context.hpp"
#include "ODEIntegrator/Context/OutFileService.hpp"
#include "ODEIntegrator/InputOutput/Reader/csvReader.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

using std::string;

void testReadGoodFile() {
    auto fileService = std::make_shared<OutFileService>();
    auto mesh = std::make_shared<Mesh>();
    auto ctx = std::make_shared<Context>(fileService, mesh);
    csvReader reader(ctx, "data/good_input_file.csv");
    reader.read();
    assert(mesh->totalSize() > 0);
}

void testReadBadFile() {
    string const expectedMsg = "Input file does not have temperature!!";
    auto fileService = std::make_shared<OutFileService>();
    auto mesh = std::make_shared<Mesh>();
    auto ctx = std::make_shared<Context>(fileService, mesh);
    csvReader reader(ctx, "data/file_without_temperature.csv");
    try {
        reader.read();
    } catch (std::runtime_error const &e) {
        assert(e.what() == expectedMsg);
    }
}

void testReadNotExistingFile() {
    string const expectedMsg = "Failing opening input file: data/nothing.csv";
    auto fileService = std::make_shared<OutFileService>();
    auto mesh = std::make_shared<Mesh>();
    auto ctx = std::make_shared<Context>(fileService, mesh);
    csvReader reader(ctx, "data/nothing.csv");
    try {
        reader.read();
    } catch (std::runtime_error const &e) {
        assert(e.what() == expectedMsg);
    }
}

void testReadTxtFile() {
    string const expectedMsg = "Input file should end with '.csv', actual: data/dummmy.txt";
    auto fileService = std::make_shared<OutFileService>();
    auto mesh = std::make_shared<Mesh>();
    auto ctx = std::make_shared<Context>(fileService, mesh);
    try {
        csvReader reader(ctx, "data/dummmy.txt");
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