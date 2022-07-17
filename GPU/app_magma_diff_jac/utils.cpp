#include <app_magma_diff_jac/utils.hpp>


namespace Utils
{

    /*
        Print number of active threads
    */
    // int printNumThreads(){
    //     int numThreads;
    //     #pragma omp parallel
    //     {
    //         int id = omp_get_thread_num();
    //         if (id == 0){
    //             numThreads = omp_get_num_threads();
    //             std::cout << "Number of threads: " << numThreads << std::endl;
    //         }
    //     }
    //     return numThreads;
    // }

    /*
        Count the number of species are in the input file

        Reference bool to indicate if coordinates data is provided:
            coordFlag = false -> No coords
            coordFlag = true -> coords need to be read
    */
    size_t numSpeciesCsv(const std::string& header, bool &coordFlag){
        size_t nsp {0};

        // Format of the species match
        std::string match {"CO"};

        // Last Column if there is no coords in the input file
        std::string noCoords {"TEM"};

        std::string item;
        std::stringstream itemIterator(header);

        // Vector of strings containing the names of the columns
        std::vector<std::string> columns;

        // Iterate the header's words and count species
        while(std::getline(itemIterator, item, ',')){

            if (item.find(match) != std::string::npos)
                nsp++;

            columns.push_back(item);
        }

        // Check if coordinates data is provided:
        coordFlag = ( columns.back().find(noCoords) == std::string::npos ) ? true : false;

        return nsp;
    }



    /*
        Read csv file with the thermodynamic data for a 
        time-step of the simulation
        
        Input: 
            - string: name of the file
        Output:
            - pointer to a struct: contains the species and therm data of each point
    */
    std::shared_ptr<ThermData> readCsv(const std::string& csvFileName){
        // Min valid mass fraction
        //double eps {1e-16};

        std::string header;

        // Read file object and open
        std::ifstream file;
        file.open(csvFileName);

        // Health check
        if (!file.is_open()){
            std::cout << "Error opening input file" << std::endl;
            exit(1);
        }


        // Get the first row because (the header)
        std::getline(file, header);
        
        // Count number of species
        bool coordFlag;
        size_t nsp = numSpeciesCsv(header, coordFlag);

        // Struct with the data of the file
        auto mesh = std::make_shared<ThermData>();
        mesh->nsp = nsp;
        mesh->header = header;
        mesh->coordFlag = coordFlag;


        /* Read dataset */
        std::string row {""};
        size_t n{0};
        while(std::getline(file, row)){
            n++;
            std::stringstream items(row);

            // Save species values in one row of the matrix
            std::vector<double> species(mesh->nsp, 0.0f);
            std::string item;
            for (double& i : species){
                std::getline(items,item,',');
                
                // Filter extremely low values of mass fraction
                // double tmp = std::stod(item);
                // if (tmp > eps)
                //     i = tmp;
                i = std::stod(item);
            }
            mesh->matSp.push_back(species);

            // Read Enthalpy and Temperature 
            std::getline(items, item, ',');
            mesh->enthal.push_back(std::stod(item));
            std::getline(items, item, ',');
            mesh->temp.push_back(std::stod(item));

            if (mesh->coordFlag){
                // Read coordinates (3D)
                Coords coords;

                std::getline(items, item, ',');
                coords.x = std::stod(item);
                std::getline(items, item, ',');
                coords.y = std::stod(item);
                std::getline(items, item, ',');
                coords.z = std::stod(item);

                mesh->pos.push_back(coords);
            }

        }

        mesh->points = n;

        // Close file
        file.close();
        return mesh;
    }


    /*
        Write csv file with the thermodynamic data after
        a time-step of the simulation
        
        Input: 
            - const shared_ptr to the struct data
            - String: name of the outputfile
    */
    void writeCsv(const std::shared_ptr<ThermData>& mesh, const std::string& csvFileName) {
        
        std::ofstream file(csvFileName);

        if (!file.is_open()){
            std::cout << "Error opening/creating output file" << std::endl;
            exit(1);
        }
        
        // Write header
        file << mesh->header << std::endl;
        
        // Write data
        size_t n = mesh->points;
        size_t nsp = mesh->nsp;
        for (size_t i = 0; i < n; i++){

            // Precision and notation
            file.precision(4);
            file << std::scientific;

            // Write species
            for (size_t j = 0; j < nsp; j++){
                file << mesh->matSp[i][j] << ',';
            }

            // Write enthalpy
            file.precision(3);
            file << mesh->enthal[i] << ',';

            // Write temperature
            file << std::defaultfloat << mesh->temp[i];

            if (mesh->coordFlag){
                // Write coordinates
                file.precision(8);
                file << ',' << mesh->pos[i].x << ',' 
                            << mesh->pos[i].y << ',' 
                            << mesh->pos[i].z;
            }

            file << std::endl;
            
        }

        file.close();
    }

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
                 const std::vector<double> time, const size_t meshPoints, bool append){
        
        std::ofstream file;

        if (!append){
            // Header of the data
            std::string header {"mechanism,threads,package,points,read_time[s],calc_time[s],write_time[s]"};

            // Overwrite/create mode
            file.open(csvFileName, std::ofstream::out);
            
            if (!file.is_open()){
                std::cout << "Error opening/creating report file" << std::endl;
                exit(1);
            }
            
            // Report
            file << header << std::endl;
            file << mechanism << ',' << numThreads << ',' << sizePackage << ',' << meshPoints;
            
            size_t n = time.size();
            for (int i = 0; i < n; i++){
                file << ',' << time[i];
            }

            file << std::endl;

        } else {
            // Append mode
            file.open(csvFileName, std::ofstream::app);

            if (!file.is_open()){
                std::cout << "Error opening/creating report file" << std::endl;
                exit(1);
            }

            // Report
            file << mechanism << ',' << numThreads << ',' << sizePackage << ',' << meshPoints;
            
            size_t n = time.size();
            for (int i = 0; i < n; i++){
                file << ',' << time[i];
            }

            file << std::endl;

        }

        file.close();
    }
}
