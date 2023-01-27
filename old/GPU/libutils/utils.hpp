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
//#include <omp.h>
#include <vector>
#include <memory>
#include <algorithm>

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
        Generic Logger
    */
    class Logger {
        public:
            Logger(std::string _header) {
                header = " " + _header + header + " ";
            }

            enum logLevels { ERROR, INFO };

            void set_logger_level(std::string _level) {
                logLevels enum_level = convert_to_log_level(_level);
                switch (enum_level)
                {
                    case INFO:
                        logLevel = INFO;
                        break;
                    
                    case ERROR:
                        logLevel = ERROR;
                        break;
                }
            }

            void print_message(std::string function, int line, std::string message, std::string log_type) {
                logLevels curr_log_type = convert_to_log_level(log_type);
                if (curr_log_type <= logLevel) {
                    std::string type = convert_log_level_to_string(curr_log_type);
                    int len = header.length();

                    std::cout << std::endl << num_char('*', 20) << header << num_char('*', 20);
                    std::cout << std::endl << std::endl;
                    std::cout << type + ": ";
                    std::cout << "From \"" + function + "\" (line: ";
                    std::cout << line;
                    std::cout << ") -> ";
                    std::cout << message << std::endl;
                    std::cout << std::endl << num_char('*', 40 + len) << std::endl << std::endl;
                }
            }


        private:
            std::string header{" logger"};
            logLevels logLevel = INFO;

            std::string num_char(char character, int times) {
                std::string chain{""};
                for(int i = 0; i < times; i++) {
                    chain.push_back(character);
                }
                return chain;
            }

            logLevels convert_to_log_level(std::string str_level) {
                // String to lowercase
                std::for_each(str_level.begin(), str_level.end(), [](char &c) {
                    c = std::tolower(c);
                });
                
                logLevels result;
                if (str_level == "info")
                    result = INFO;
                if (str_level == "error")
                    result = ERROR;

                return result;
            }

            std::string convert_log_level_to_string(logLevels level) {
                std::string result;
                if (level == INFO)
                    result = "INFO";
                if (level == ERROR)
                    result = "ERROR";
                return result;
            }
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
    //int printNumThreads();


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

#endif


