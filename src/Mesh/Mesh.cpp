#include "Mesh.hpp"

#include <vector>
#include <stdexcept>

using std::vector;

Mesh& Mesh::get() {
    static Mesh mesh;
    return mesh;
}

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

void Mesh::addPoint(Point const& newPoint) {
    if (totalSize() == 0) {
        enthalpyFlag = newPoint.hasEnthalpy();
        coordFlag = newPoint.hasCoordinates();
        nsp = newPoint.numSpecies();
    }
    if (isCompatible(newPoint)) {
        temperature.push_back(newPoint.getTemperature());
        // Move syntax to avoid copying
        speciesMatrix.push_back(std::move(newPoint.getSpecies()));

        if (hasEnthalpy())
            enthalpy.push_back(newPoint.getEnthalpy());
        if (hasCoordinates())
            coordinates.push_back(newPoint.getCoordinates());

    } else {
        throw std::runtime_error("Trying to add an incompatible Point to Mesh");
    }
}

int Mesh::totalSize() const {
    return temperature.size();
}

int Mesh::numSpecies() const {
    return nsp;
}

bool Mesh::hasEnthalpy() const {
    return enthalpyFlag;
}

bool Mesh::hasCoordinates() const {
    return coordFlag;
}

void Mesh::clear() {
    temperature.clear();
    enthalpy.clear();
    speciesMatrix.clear();
    coordinates.clear();
    enthalpyFlag = false;
    coordFlag = false;
    nsp = 0;
}

bool Mesh::isCompatible(Point const& point) const {
    return point.hasEnthalpy() == hasEnthalpy() &&
            point.hasCoordinates() == hasCoordinates() &&
            point.numSpecies() == numSpecies();
}