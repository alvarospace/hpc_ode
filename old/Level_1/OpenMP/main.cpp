// This program solves equations of type M y'= f(y) 

#include "cantera/thermo.h"
#include "cantera/zerodim.h"
#include "utils.hpp"
//using namespace Cantera;

void printNumThreads();

void run(const std::string& mechanism, const std::string& phase,
         const std::string& inputFile, const std::string& outputFile,
         const size_t pack_size, const bool append){

    /********** INPUT CONSTANTS ************/

    // Input Constants (Pa and K)
    const double p = 101325.15;

    // Time step
    const double dt = 1e-3;

    /****************************/

    // TIMERS
    Utils::Timer readTime, calculationTime, writeTime;


    /*********** READ CSV FILE and SET PACKAGE SIZE ************/
    readTime.tic();
    std::shared_ptr<Utils::ThermData> mesh = Utils::readCsv(inputFile);
    readTime.toc();
    mesh->info();
    
    // Number of points of the mesh
    size_t n_size = mesh->points;

    // Outer loop iterations
    size_t mod = n_size % pack_size;
    size_t ext_it = n_size / pack_size;
    // If the division is not exact, add 1 extra iteration (the remainder)
    if (mod > 0)
        ext_it++;

    std::cout << "Package size: " << pack_size << std::endl;
    std::cout << "Number of packages: " << ext_it << std::endl;
    std::cout << "Remainder points (will be process by the last package): " << mod << std::endl;
    /************************************************************/

    int nThreads = Utils::printNumThreads();


    calculationTime.tic();
    /******************** CALCULATION **********************/

    /* Outer loop to deliver work to threads */
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < ext_it; i++){
        size_t thread_points;
        // Limit the pack_size to the last iteration
        // because the number of points to process is the remainder
        if (i == (ext_it - 1) && mod > 0)
            thread_points = mod;
        else
            thread_points = pack_size;

        /********* VECTORS ***********/
        // Allocation of the required memory for Cantera application

        // Vector of pointers for solutions, ThermoPhases and reactors
        std::vector<std::shared_ptr<Cantera::Solution>> solutions(thread_points, Cantera::newSolution(mechanism,phase,"None"));
        std::vector<std::shared_ptr<Cantera::ThermoPhase>> gases(thread_points);
        std::vector<Cantera::IdealGasConstPressureReactor> reactors(thread_points);
        std::vector<Cantera::ReactorNet> nets(thread_points);
        /****************************/
        

        /* Inner loop to integrate the chemical */
        for (size_t j = 0; j < thread_points; j++){
            // Index of the global point
            size_t index = i * pack_size + j;

            // Define initial Values of unknows (Corresponds to phi = 1.0, perfect fuel air mixture)
            gases[j] = solutions[j]->thermo();
            gases[j]->setState_TPY(mesh->temp[index], p, mesh->matSp[index].data());


            // Insert the reaction
            reactors[j].insert(solutions[j]);
            nets[j].addReactor(reactors[j]);

            // Integrate time step
            nets[j].advance(dt);    

            // Save solution
            mesh->temp[index] = reactors[j].temperature();
            gases[j]->getMassFractions(mesh->matSp[index].data());
            mesh->enthal[index] = gases[j]->enthalpy_mass();
        }
    }
    /************************************************************/
    calculationTime.toc();


    writeTime.tic();
    /***** Write Solution *******/

    Utils::writeCsv(mesh, outputFile);

    /*****************************/
    writeTime.toc();


    /* Time report */
    short w1{15}, w2{20};
    std::cout.precision(8);
    std::cout << std::endl << std::fixed;
    std::cout << "\t*** TIME REPORT ***" << std::endl << std::endl;
    std::cout << std::setw(w1) << "Read:" << std::setw(w2) << readTime.time() << " s" << std::endl;
    //std::cout << std::setw(w1) << "Allocation:" << std::setw(w2) << allocateTime.time() << " ms" << std::endl;
    std::cout << std::setw(w1) << "Calculation:" << std::setw(w2) << calculationTime.time() << " s" << std::endl;
    std::cout << std::setw(w1) << "Write:" << std::setw(w2) << writeTime.time() << " s" << std::endl;
    double totalTime = readTime.time() + calculationTime.time() + writeTime.time();
    std::cout << std::endl << std::setw(w1) << "TOTAL TIME:" << std::setw(w2) << totalTime << " s" << std::endl;

    /* Write csv time report */
    std::vector<double> time = {readTime.time(), calculationTime.time(), writeTime.time()};
    Utils::reportCsv("report_time.csv", mechanism, nThreads, pack_size, time, n_size, append);

}


int main(int argc, char *argv[]){
    std::string mechanism {"gri30.yaml"};
    std::string phase {""};
    std::string inputFile {"res2000.csv"};
    std::string outputFile {"out.csv"};
    size_t pack{10};
    bool append{false};

    /* Command-line arguments logic */
    std::vector<std::string> args;
    if (argc > 1) {
        for (int i = 0; i < argc; i++) {
            args.push_back(std::string(argv[i]));
        }

        if (args[1].compare("--help") == 0 || args[1].compare("-h") == 0) {
            std::cout << "Usage: " << args[0] <<" [input file=res2000.csv] [output file=out.csv] [mechanism=gri30.yaml] [package size=10] [append report file=0]" << std::endl;
            return 0;
        } else {
            if (argc == 2) {
                inputFile = args[1];
            } else if (argc == 3) {
                inputFile = args[1];
                outputFile = args[2];
            } else if (argc == 4) {
                inputFile = args[1];
                outputFile = args[2];
                mechanism = args[3];
            } else if (argc == 5) {
                inputFile = args[1];
                outputFile = args[2];
                mechanism = args[3];
                pack = std::stoi(args[4]);
            } else if (argc == 6) {
                inputFile = args[1];
                outputFile = args[2];
                mechanism = args[3];
                pack = std::stoi(args[4]);
                if (std::stoi(args[5]) == 1){
                    append = true;
                }
            } else {
                std::cout << "Too many arguments..." << std::endl;
                return 1;
            }
        }
    }

    /* Program */
    try {
        run(mechanism, phase, inputFile, outputFile, pack, append);
    } catch (Cantera::CanteraError& err) {
        std::cout << err.what() << std::endl;
        Cantera::appdelete();
        return 1;
    }
    return 0;
}


