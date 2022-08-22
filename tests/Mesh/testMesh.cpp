#include "Mesh/Mesh.hpp"
#include "Mesh/Point.hpp"
#include "Mesh/PointBuilder.hpp"
#include <cassert>
#include <array>
#include <stdexcept>

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

Point createBadPoint() {
    int const nsp {3};
    vector<double> species {0.2, 0.2, 0.6};

    Point point = Point::create(nsp)
                        .addSpecies(species);

    return point;
}

void testAddPointAndClearMesh() {
    // Add 3 points to the mesh
    array<int,3> range {0, 1, 2};

    Mesh &mesh = Mesh::get();
    for (auto const& i : range) {
        Point point = createPoint();
        mesh.addPoint(point);
    }

    assert(mesh.totalSize() == 3);
    mesh.clear();
    assert(mesh.totalSize() == 0);
}

void testOutOfIndexAccessToMesh() {
    // Add 3 points to the mesh
    array<int,3> range {0, 1, 2};

    Mesh &mesh = Mesh::get();
    for (auto const& i : range) {
        Point point = createPoint();
        mesh.addPoint(point);
    }

    vector<double> v1 = mesh.getSpeciesVector(1);
    string const expectedMsg = "vector::_M_range_check: __n (which is 4) >= this->size() (which is 3)";
    try {
        vector<double> v4 = mesh.getSpeciesVector(4);
    } catch (out_of_range const &e) {
        string const errorMsg = e.what();
        assert(errorMsg == expectedMsg);
    }
    
    mesh.clear();
}

void testAddBadPointToMesh() {
    string const expectedMsg = "Trying to add an incompatible Point to Mesh";
    Mesh &mesh = Mesh::get();
    mesh.clear();
    Point badPoint = createBadPoint();
    try {
        mesh.addPoint(badPoint);
    } catch (exception const &e) {
        string const errorMsg = e.what();
        assert(errorMsg == expectedMsg);
    }
    assert(mesh.totalSize() == 0);
}

void AddPointToMesh() {
    Mesh &mesh = Mesh::get();
    Point point = createPoint();
    mesh.addPoint(point);
}

void testMeshState() {
    AddPointToMesh();
    Mesh &mesh = Mesh::get();
    Point point = createPoint();
    mesh.addPoint(point);
    assert(mesh.totalSize() == 2);
    mesh.clear();
}

int main() {
    testAddPointAndClearMesh();
    testOutOfIndexAccessToMesh();
    testAddBadPointToMesh();
    testMeshState();
    return 0;
}