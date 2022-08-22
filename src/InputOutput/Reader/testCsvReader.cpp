#include "csvReader.hpp"

void testRead() {
    csvReader reader("hola.txt");
    reader.read();
}

int main() {
    return 0;
}