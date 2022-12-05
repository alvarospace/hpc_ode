#pragma once

#include <vector>

#include "ODEIntegrator/Mesh/Point.hpp"

// Class that holds the data of the execution
class Mesh {
    public:
        Mesh() {}
        
        std::vector<double>& getTemperatureVector();
        std::vector<double>& getEnthalpyVector();
        std::vector<double>& getSpeciesVector(int index);
        std::vector<std::vector<double>>& getSpeciesMatrix();
        std::vector<Coords>& getCoordinatesVector();

        double* getTemperaturePointer();
        double* getEnthalpyPointer();
        double* getSpeciesPointer(int index);

        void addPoint(Point const& newPoint);
        int totalSize() const;
        int numSpecies() const;
        bool hasEnthalpy() const;
        bool hasCoordinates() const;
        void clear();

        // Remove copy constructor and assignment operator
        Mesh(Mesh const&) = delete;
        Mesh& operator=(Mesh const&) = delete;
    
    private:
        bool isCompatible(Point const& point) const;

        std::vector<double> temperature;
        std::vector<double> enthalpy;
        std::vector<std::vector<double>> speciesMatrix;
        std::vector<Coords> coordinates;

        int nsp {0};
        bool enthalpyFlag {false};
        bool coordFlag {false};
};