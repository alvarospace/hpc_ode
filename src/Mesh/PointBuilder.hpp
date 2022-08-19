#pragma once

#include "Point.hpp"
#include <vector>

using std::vector;

class PointBuilder {
    public:
        PointBuilder(int _nsp) : point(Point(_nsp)) {}

        // Needed operator for the builder pattern
        operator Point() const { return std::move(point); }

        // Build steps
        PointBuilder& addTemperature(double temperature);
        PointBuilder& addEnthalpy(double enthalpy);
        PointBuilder& addSpecies(vector<double> species);
        PointBuilder& addCoordinates(Coords coords);


    private:
        Point point;
};