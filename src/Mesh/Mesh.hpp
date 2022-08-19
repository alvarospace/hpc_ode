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

        void addPoint();
        int size() const;
        bool hasCoords() const;
    
    private:
        vector<double> temperature;
        vector<double> enthalpy;
        vector<vector<double>> speciesMatrix;
        vector<Coords> coordinates;

        int Points;
        int nsp;
        bool coordFlag;
};