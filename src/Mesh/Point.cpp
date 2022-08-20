#include "Point.hpp"
#include "PointBuilder.hpp"
#include <stdexcept>
#include <vector>

using std::vector;

PointBuilder Point::create(int _nsp) {
    return PointBuilder{_nsp};
}

bool Point::hasTemperature() const {
    return flagTemperature;
}

bool Point::hasEnthalpy() const {
    return flagEnthapy;
}

bool Point::hasSpecies() const {
    return flagSpecies;
}

bool Point::hasCoordinates() const {
    return flagCoordinates;
}

bool Point::isReady() const {
    return flagTemperature && flagSpecies;
}

double Point::getTemperature() const {
    if (!hasTemperature()) {
        throw std::logic_error("This Point does not have temperature");
    }
    return temperature;
}

double Point::getEnthalpy() const {
    if (!hasEnthalpy()) {
        throw std::logic_error("This Point does not have enthalpy");
    }
    return enthalpy;
}

vector<double> Point::getSpecies() const {
    if (!hasSpecies()) {
        throw std::logic_error("This Point does not have species");
    }
    return species;
}

Coords Point::getCoordinates() const {
    if (!hasCoordinates()) {
        throw std::logic_error("This Point does not have coordinates");
    }
    return coordinates;
}