#include "ODEIntegrator/Mesh/Mesh.hpp"

#include <vector>
#include <stdexcept>

using std::vector;

vector<double>& Mesh::getTemperatureVector() {
    return temperature;
}

vector<double>& Mesh::getEnthalpyVector() {
    if (!hasEnthalpy())
        throw std::runtime_error("Mesh does not have enthalpy");
    return enthalpy;
}

vector<double>& Mesh::getSpeciesVector(int index) {
    return speciesMatrix.at(index);
}

vector<vector<double>>& Mesh::getSpeciesMatrix() {
    return speciesMatrix;
}

vector<Coords>& Mesh::getCoordinatesVector() {
    if (!hasCoordinates())
        throw std::runtime_error("Mesh does not have coordinates");
    return coordinates;
}

double* Mesh::getTemperaturePointer() {
    return temperature.data();
}

double* Mesh::getEnthalpyPointer() {
    if (!hasEnthalpy())
        throw std::runtime_error("Mesh does not have enthalpy");
    return enthalpy.data();
}

double* Mesh::getSpeciesPointer(int index) {
    return speciesMatrix.at(index).data();
}

void Mesh::addPoint(Point const& newPoint) {
    if (totalSize() == 0 && newPoint.isReady()) {
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
            point.numSpecies() == numSpecies() &&
            point.isReady();
}