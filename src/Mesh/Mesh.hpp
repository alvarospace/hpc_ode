#pragma once
#include "Point.hpp"
#include <vector>

using std::vector;

// Has to be a Singleton
class Mesh {
    public:
        vector<double>& getTemperatureVector();
        vector<double>& getEnthalpyVector();
        vector<double>& getSpeciesVector(int index);
        vector<vector<double>>& getSpeciesMatrix();
        vector<Coords>& getCoordinatesVector();

        void addPoint(Point const& newPoint);
        int size() const;
        bool hasEnthalpy() const;
        bool hasCoordinates() const;
    
    private:
        bool isCompatible(Point const& point) const;

        vector<double> temperature;
        vector<double> enthalpy;
        vector<vector<double>> speciesMatrix;
        vector<Coords> coordinates;

        int nsp;
        bool enthalpyFlag {false};
        bool coordFlag {false};
};