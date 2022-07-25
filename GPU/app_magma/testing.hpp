#ifndef TEST

#define TEST
#include <memory>
#include <iostream>
#include <string>
#include <map>

#include <cuda.h>
#include <cuda_runtime.h>
#include <libutils/utils.hpp>

using namespace std;

namespace Testing
{
    class YsunYpyjac {

        public:
            YsunYpyjac(double *_ySunHost, double *_ySunGPU, double *_yPyjacHost,
                        double *_yPyjacGPU, double *_dyPyjacHost, double *_dyPyjacGPU);

            void set_simulation_size(int num_systems, int num_species, int _padded);

            void set_logger_level(string level);

            // Interface 1 of compare_arrays
            void ysun_vs_ypyjac();

            // Interface 2 of compare_arrays
            void ysun_vs_dypyjac();

        private:
            int    systems;
            int    nsp;
            int    padded;
            int    total_bytes_sun;
            int    total_bytes_py;
            double *ySunHost {nullptr};
            double *ySunGPU {nullptr};
            double *yPyjacHost {nullptr};
            double *yPyjacGPU {nullptr};
            double *dyPyjacHost {nullptr};
            double *dyPyjacGPU {nullptr};
            unique_ptr<Utils::Logger> logger;
            const string info_type = "info";
            const string error_type = "error";
            map<string, double*> pointer_dict;

            void compare_arrays(string array1Host, string array1GPU, string array2Host, string array2GPU);
    };
}

#endif