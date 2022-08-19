#pragma once
#include "System.h"
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

        void addSystem();
        int size() const;
        bool hasCoords() const;
    
    private:
        vector<double> temperature;
        vector<double> enthalpy;
        vector<vector<double>> speciesMatrix;
        vector<Coords> coordinates;

        int systems;
        int nsp;
        bool coordFlag;
};