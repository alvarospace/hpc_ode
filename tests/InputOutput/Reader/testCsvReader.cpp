#include "Reader/csvReader.hpp"

void testRead() {
    csvReader reader("data/csv_example.csv");
    reader.read();
}

int main() {
    testRead();
    return 0;
}