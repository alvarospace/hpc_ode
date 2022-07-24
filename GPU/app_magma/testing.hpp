#ifndef TEST

#define TEST
#include <memory>
#include <iostream>
#include <string>

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

            // TODO: abstract layer to evaluate both "y vs y" or "y vs dy" 
            void ysun_vs_ypyjac();

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
            unique_ptr<Utils::Logger> infoLogger;
            unique_ptr<Utils::Logger> errorLogger;
            string logLevel{"error"};
    };
}

#endif