#pragma once

#include "System.h"
#include <vector>

using std::vector;

class SystemBuilder {
    public:
        SystemBuilder(int _nsp) : system(System(_nsp)) {}

        // Needed operator for the builder pattern
        operator System() const { return std::move(system); }

        // Build steps
        SystemBuilder& addTemperature(double temperature);
        SystemBuilder& addEnthalpy(double enthalpy);
        SystemBuilder& addSpecies(vector<double> species);
        SystemBuilder& addCoordinates(Coords coords);


    private:
        System system;
};