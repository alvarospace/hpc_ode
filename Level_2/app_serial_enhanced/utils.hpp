/*
    Small library with utilities for TFM-BSC project
*/

#ifndef UTILS

#define UTILS
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <vector>
#include <memory>

#endif

namespace Utils
{   
    /*
        User defined class to measure time
    */
    class Timer {
        public:
            void tic() {
                t1 = std::chrono::steady_clock::now();
            }

            void toc() {
                t2 = std::chrono::steady_clock::now();
                measured_time = t2 - t1;
            }

            double time() {
                return measured_time.count();
            }


        private:
            std::chrono::steady_clock::time_point t1;
            std::chrono::steady_clock::time_point t2;
            std::chrono::duration<double> measured_time;
    };




    /*
        Struct to save the spacial coordinates of a point in a mesh
    */
    struct Coords {
        double x,y,z;
    };

    /*
        Struct that contains the thermodynamic data of each point of the input file
    */
    struct ThermData {
        size_t points{0}, nsp{0};
        std::vector<double> temp;
        std::vector<double> enthal;
        std::vector<std::vector<double>> matSp;
        std::vector<Coords> pos;
        std::string header{""};
        bool coordFlag;

        void info(){
            std::cout << std::endl << "********* ThermData INFO **********" << std::endl << std::endl;

            std::cout << "Number of species found in the csv file: " << nsp << std::endl;
            std::cout << "Number of points: " << points << std::endl;

            if (!coordFlag){
                std::cout << "Coords points NOT FOUND" << std::endl;
            }

            std::cout << std::endl << "***********************************" << std::endl << std::endl;
        }
    };


    /*
        Print number of active threads
    */
    int printNumThreads();


    /*
        Count the number of species are in the input file

        Reference bool to indicate if coordinates data is provided:
            coordFlag = false -> No coords
            coordFlag = true -> coords need to be read
    */
    size_t numSpeciesCsv(const std::string& header, bool &coordFlag);


    /*
        Read csv file with the thermodynamic data for a 
        time-step of the simulation
        
        Input: 
            - string: name of the file
        Output:
            - pointer to a struct: contains the species and therm data of each point
    */
    std::shared_ptr<ThermData> readCsv(const std::string& csvFileName);


    /*
        Write csv file with the thermodynamic data after
        a time-step of the simulation
        
        Input: 
            - const shared_ptr to the struct data
            - String: name of the outputfile
    */
    void writeCsv(const std::shared_ptr<ThermData>& mesh, const std::string& csvFileName);


    /*
        Write Csv file with the time report of the execution, the objetive is the analysis
        the number of OpenMP threads and package size with best performance
        
        Input:
            - filename
            - number of threads
            - size of the package processed by thread
            - vector of double containing the time report
            - number of points in the input file
            - boolean to append report or create new file
    */
    void reportCsv(const std::string& csvFileName, const std::string& mechanism, const int numThreads, const size_t sizePackage,
                 const std::vector<double> time, const size_t meshPoints, bool append);


}
