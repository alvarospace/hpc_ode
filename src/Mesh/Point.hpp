#pragma once

#include "PointBuilder.hpp"
#include <vector>

using std::vector;

struct Coords {
    double x,y,z;
};

// Implements builder pattern
class Point {
    public:
        Point(int _nsp) : nsp(_nsp) {}

        friend class PointBuilder;
        static PointBuilder create(int _nsp);

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