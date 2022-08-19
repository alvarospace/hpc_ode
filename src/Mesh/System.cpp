#include "System.h"
#include "SystemBuilder.h"
#include <stdexcept>
#include <vector>

using std::vector;

SystemBuilder System::create(int _nsp) {
    return SystemBuilder{_nsp};
}

bool System::hasTemperature() {
    return flagTemperature;
}

bool System::hasEnthalpy() {
    return flagEnthapy;
}

bool System::hasSpecies() {
    return flagSpecies;
}

bool System::hasCoordinates() {
    return flagEnthapy;
}

double System::getTemperature() {
    if (!hasTemperature()) {
        throw std::logic_error("This system does not have temperature");
    }
    return temperature;
}

double System::getEnthalpy() {
    if (!hasEnthalpy()) {
        throw std::logic_error("This system does not have enthalpy");
    }
    return enthalpy;
}

vector<double> System::getSpecies() {
    if (!hasSpecies()) {
        throw std::logic_error("This system does not have species");
    }
    return species;
}

Coords System::getCoordinates() {
    if (!hasCoordinates) {
        throw std::logic_error("This system does not have coordinates");
    }
    return coordinates;
}