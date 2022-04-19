#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>

using std::ofstream;
using namespace std;

int main()
{
   string fname="res.csv";
   string line, word;
   double data;
   fstream file (fname, ios::in);
   while(getline(file, line)) {
        row.clear();
        stringstream str(line);
        while(getline(str, word, ',')){
             data = stod(word);
             row.push_back(data);
	}
        content.push_back(row);
   }
}
