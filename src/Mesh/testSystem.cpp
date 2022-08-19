#include "Point.hpp"
#include "PointBuilder.hpp"
#include <iostream>

using namespace std;

// Dummy values for testing
int const nsp {3};
double const temperature {100};
double const enthalpy {-10000};
Coords const coordinates {1,2,3};
vector<double> species {0.2, 0.2, 0.6};

void testFullBuild() {
    
    Point point = Point::create(nsp)
                            .addTemperature(temperature)
                            .addEnthalpy(enthalpy)
                            .addCoordinates(coordinates)
                            .addSpecies(species);
}





int main() {
    testFullBuild();
    return 0;
}