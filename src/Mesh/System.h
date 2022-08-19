#pragma once

#include "SystemBuilder.h"
#include <vector>

using std::vector;

struct Coords {
    double x,y,z;
};

// Implements builder pattern
class System {
    public:
        System(int _nsp) : nsp(_nsp) {}

        friend class SystemBuilder;
        static SystemBuilder create(int _nsp);

        bool hasTemperature();
        bool hasEnthalpy();
        bool hasSpecies();
        bool hasCoordinates();

        double getTemperature();
        double getEnthalpy();
        vector<double> getSpecies();
        Coords getCoordinates();

    private:
        double temperature;
        double enthalpy;
        vector<double> species;
        Coords coordinates;
        int nsp;

        bool flagTemperature{false};
        bool flagEnthapy{false};
        bool flagSpecies{false};
        bool flagCoordinates{false};
};