#include "System.h"
#include <iostream>

using namespace std;

void testFullBuild() {
    int const nsp {3};
    double const temperature {100};
    double const enthalpy {-10000};
    Coords const coordinates {1,2,3};

    cout << "x: " << coordinates.x << endl;
    cout << "y: " << coordinates.y << endl;
    cout << "z: " << coordinates.z << endl;
}





int main() {
    testFullBuild();
    return 0;
}