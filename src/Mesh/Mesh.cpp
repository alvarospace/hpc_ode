#include "Mesh.hpp"

#include <vector>
#include <stdexcept>

using std::vector;

vector<double>& Mesh::getTemperatureVector() {
    return temperature;
}

vector<double>& Mesh::getEnthalpyVector() {
    return enthalpy;
}

vector<double>& Mesh::getSpeciesVector(int index) {
    return speciesMatrix.at(index);
}

vector<vector<double>>& Mesh::getSpeciesMatrix() {
    return speciesMatrix;
}

vector<Coords>& Mesh::getCoordinatesVector() {
    return coordinates;
}