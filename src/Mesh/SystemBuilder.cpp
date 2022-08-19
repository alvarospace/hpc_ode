#include "SystemBuilder.h"
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <numeric>

SystemBuilder& SystemBuilder::addTemperature(double temperature) {
    if (system.hasTemperature()) {
        throw std::logic_error("Temperature already set in the build process");
    }
    if (temperature < 0) {
        throw std::runtime_error("Temperature cannot be negative");
    }
    system.temperature = temperature;
    system.flagTemperature = true;
    return *this;
}

SystemBuilder& SystemBuilder::addEnthalpy(double enthalpy) {
    if (system.hasEnthalpy()) {
        throw std::logic_error("Enthalpy already set in the build process");
    }
    system.enthalpy = enthalpy;
    system.flagEnthapy = true;
    return *this;
}

SystemBuilder& SystemBuilder::addSpecies(vector<double> species) {
    if (system.hasSpecies()) {
        throw std::logic_error("Species vector already set in the build process");
    }
    if (species.size() != system.nsp) {
        throw std::runtime_error("Vector size does not match the "
                "number of the species of the system");
    }

    bool all_less_equal_1 = std::all_of(std::begin(species), std::end(species), [] (double sp) {
        return sp <= 1.0;
    });
    
    int total_sum = std::accumulate(std::begin(species), std::end(species), 0.0);
    double error = 0.01;
    bool mass_conservation = total_sum > (1.0 - error) && total_sum < (1.0 + error);

    if (!all_less_equal_1 || !mass_conservation) {
        throw std::runtime_error("Species vector breaks mass conservation law");
    }

    system.species = species;
    system.flagSpecies = true;
    return *this;
}

SystemBuilder& SystemBuilder::addCoordinates(Coords coordinates) {
    if (system.hasCoordinates()) {
        throw std::logic_error("Coordinates already set in the build process");
    }
    system.coordinates = coordinates;
    system.flagCoordinates = true;
    return *this;
}