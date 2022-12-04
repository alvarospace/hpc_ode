#include <string>
#include <cassert>

#include "ODEIntegrator/Mesh/Point.hpp"
#include "ODEIntegrator/Mesh/PointBuilder.hpp"

using namespace std;

// Dummy values for testing
int const nsp {3};
double const temperature {100};
double const enthalpy {-10000};
Coords const coordinates {1,2,3};
vector<double> species {0.2, 0.2, 0.6};
vector<double> badSpecies {1.2, 0.5, 7};

void testFullBuild() {
    Point point = Point::create(nsp)
                        .addTemperature(temperature)
                        .addEnthalpy(enthalpy)
                        .addCoordinates(coordinates)
                        .addSpecies(species);
    assert(point.hasTemperature());
    assert(point.hasEnthalpy());
    assert(point.hasCoordinates());
    assert(point.hasSpecies());
}

void testMinBuild() {
    Point point = Point::create(nsp)
                        .addTemperature(temperature)
                        .addSpecies(species);
    assert(point.isReady());
}

void testBadBuild() {
    Point point = Point::create(nsp)
                        .addTemperature(temperature)
                        .addEnthalpy(enthalpy);
    assert(!point.isReady());
}

void testBadConstruction() {
    string const expected_msg = "Temperature already set in the build process";
    string error_msg;
    try {
        Point point = Point::create(nsp)
                            .addTemperature(temperature)
                            .addEnthalpy(enthalpy)
                            .addCoordinates(coordinates)
                            .addSpecies(species)
                            .addTemperature(temperature);
    } catch (exception &e) {
        error_msg = e.what();
    }
    assert(error_msg == expected_msg);
}

void testBadSpecies() {
    string const expected_msg = "Species vector breaks mass conservation law";
    string error_msg;
    try {
        Point point = Point::create(nsp)
                            .addTemperature(temperature)
                            .addEnthalpy(enthalpy)
                            .addCoordinates(coordinates)
                            .addSpecies(badSpecies);
    } catch (exception &e) {
        error_msg = e.what();
    }
    assert(error_msg == expected_msg);
}

int main() {
    testFullBuild();
    testMinBuild();
    testBadBuild();
    testBadConstruction();
    testBadSpecies();
    return 0;
}