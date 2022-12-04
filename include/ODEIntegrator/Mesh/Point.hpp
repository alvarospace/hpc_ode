#pragma once

#include <vector>

struct Coords {
    double x,y,z;
};

// Forward declaration to avoid import conflicts
class PointBuilder;

// Implements builder pattern
class Point {
    public:
        friend class PointBuilder;
        static PointBuilder create(int _nsp);

        bool hasTemperature() const;
        bool hasEnthalpy() const;
        bool hasSpecies() const;
        bool hasCoordinates() const;
        bool isReady() const;

        double getTemperature() const;
        double getEnthalpy() const;
        std::vector<double> getSpecies() const;
        Coords getCoordinates() const;
        int numSpecies() const;

    private:
        // Private constructor to force PointBuilder functionality
        Point(int _nsp) : nsp(_nsp) {}

        double temperature;
        double enthalpy;
        std::vector<double> species;
        Coords coordinates;
        int nsp;

        bool flagTemperature{false};
        bool flagEnthapy{false};
        bool flagSpecies{false};
        bool flagCoordinates{false};
};