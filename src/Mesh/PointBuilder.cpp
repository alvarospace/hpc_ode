#include "Mesh/PointBuilder.hpp"
#include "Mesh/Point.hpp"
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <numeric>

PointBuilder& PointBuilder::addTemperature(double temperature) {
    if (point.hasTemperature()) {
        throw std::logic_error("Temperature already set in the build process");
    }
    if (temperature < 0) {
        throw std::runtime_error("Temperature cannot be negative");
    }
    point.temperature = temperature;
    point.flagTemperature = true;
    return *this;
}

PointBuilder& PointBuilder::addEnthalpy(double enthalpy) {
    if (point.hasEnthalpy()) {
        throw std::logic_error("Enthalpy already set in the build process");
    }
    point.enthalpy = enthalpy;
    point.flagEnthapy = true;
    return *this;
}

PointBuilder& PointBuilder::addSpecies(vector<double> species) {
    if (point.hasSpecies()) {
        throw std::logic_error("Species vector already set in the build process");
    }
    if (species.size() != point.numSpecies()) {
        throw std::runtime_error("Vector size does not match the "
                "number of the species of the Point");
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

    point.species = species;
    point.flagSpecies = true;
    return *this;
}

PointBuilder& PointBuilder::addCoordinates(Coords coordinates) {
    if (point.hasCoordinates()) {
        throw std::logic_error("Coordinates already set in the build process");
    }
    point.coordinates = coordinates;
    point.flagCoordinates = true;
    return *this;
}