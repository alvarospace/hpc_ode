#pragma once
#include "ODEIntegrator/Mesh/Point.hpp"
#include <vector>

using std::vector;

// Singleton class to keep the application data
class Mesh {
    public:
        // Static method with static instace of Mesh
        static Mesh& get();

        vector<double>& getTemperatureVector();
        vector<double>& getEnthalpyVector();
        vector<double>& getSpeciesVector(int index);
        vector<vector<double>>& getSpeciesMatrix();
        vector<Coords>& getCoordinatesVector();

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
        // Singleton constructor is private
        Mesh() {}
        bool isCompatible(Point const& point) const;

        vector<double> temperature;
        vector<double> enthalpy;
        vector<vector<double>> speciesMatrix;
        vector<Coords> coordinates;

        int nsp {0};
        bool enthalpyFlag {false};
        bool coordFlag {false};
};