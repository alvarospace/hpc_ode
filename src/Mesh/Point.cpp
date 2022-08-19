#include "Point.hpp"
#include "PointBuilder.hpp"
#include <stdexcept>
#include <vector>

using std::vector;

PointBuilder Point::create(int _nsp) {
    return PointBuilder{_nsp};
}

bool Point::hasTemperature() {
    return flagTemperature;
}

bool Point::hasEnthalpy() {
    return flagEnthapy;
}

bool Point::hasSpecies() {
    return flagSpecies;
}

bool Point::hasCoordinates() {
    return flagCoordinates;
}

double Point::getTemperature() {
    if (!hasTemperature()) {
        throw std::logic_error("This Point does not have temperature");
    }
    return temperature;
}

double Point::getEnthalpy() {
    if (!hasEnthalpy()) {
        throw std::logic_error("This Point does not have enthalpy");
    }
    return enthalpy;
}

vector<double> Point::getSpecies() {
    if (!hasSpecies()) {
        throw std::logic_error("This Point does not have species");
    }
    return species;
}

Coords Point::getCoordinates() {
    if (!hasCoordinates()) {
        throw std::logic_error("This Point does not have coordinates");
    }
    return coordinates;
}