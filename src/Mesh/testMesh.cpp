#include "Mesh.hpp"
#include "Point.hpp"
#include "PointBuilder.hpp"
#include <cassert>
#include <iostream>

using namespace std;

// Same point always (Setup)
Point createPoint() {
    int const nsp {3};
    double const temperature {100};
    double const enthalpy {-10000};
    Coords const coordinates {1,2,3};
    vector<double> species {0.2, 0.2, 0.6};

    Point point = Point::create(nsp)
                        .addTemperature(temperature)
                        .addSpecies(species)
                        .addEnthalpy(enthalpy)
                        .addCoordinates(coordinates);

    return point;
}

void testAddPointToMesh() {
    // Add 3 points to the mesh
    int const numPoints {3};
    vector<int> range(numPoints, 0);
    Mesh &mesh = Mesh::get();
    for (auto const& i : range) {
        Point point = createPoint();
        mesh.addPoint(point);
    }

}


int main() {
    testAddPointToMesh();
    return 0;
}