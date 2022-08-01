#include <app_magma/testing.hpp>

using namespace std;

namespace Testing
{
    /************************* YsunYpyjac ***************************/
    /* PUBLIC */

    // Constructor
    YsunYpyjac::YsunYpyjac(double *_ySunHost, double *_ySunGPU, double *_yPyjacHost, double *_yPyjacGPU, double *_dyPyjacHost, double *_dyPyjacGPU) \
        : ySunHost(_ySunHost), ySunGPU(_ySunGPU), yPyjacHost(_yPyjacHost), yPyjacGPU(_yPyjacGPU), dyPyjacHost(_dyPyjacHost), dyPyjacGPU(_dyPyjacGPU) {

        logger = make_unique<Utils::Logger>("Y array");

        // Dictionary of pointers
        pointer_dict["ySunHost"]     = ySunHost;
        pointer_dict["ySunGPU"]      = ySunGPU;
        pointer_dict["yPyjacHost"]   = yPyjacHost;
        pointer_dict["yPyjacGPU"]    = yPyjacGPU;
        pointer_dict["dyPyjacHost"]  = dyPyjacHost;
        pointer_dict["dyPyjacGPU"]   = dyPyjacGPU;
    }

    void YsunYpyjac::set_simulation_size(int num_systems, int num_species, int _padded) {
        systems = num_systems;
        nsp = num_species;

        // Padded is the number of systems used by pyjac
        padded = _padded;

        total_bytes_sun = systems * nsp * sizeof(double);
        total_bytes_py = padded * nsp * sizeof(double);
    }

    /* Level can be "error" or "info" */
    void YsunYpyjac::set_logger_level(string level) {
        logger->set_logger_level(level);
    }

    void YsunYpyjac::ysun_vs_ypyjac() {
        compare_arrays("ySunHost", "ySunGPU", "yPyjacHost", "yPyjacGPU");
    }

    void YsunYpyjac::ysun_vs_dypyjac() {
        compare_arrays("ySunHost", "ySunGPU", "dyPyjacHost", "dyPyjacGPU");
    }

    /* PRIVATE */

    void YsunYpyjac::compare_arrays(string array1Host, string array1GPU, string array2Host, string array2GPU) {
        // Find pointers in dictionary
        double *ptr_array1Host;
        double *ptr_array1GPU;
        double *ptr_array2Host;
        double *ptr_array2GPU;
        try {
            ptr_array1Host = pointer_dict.at(array1Host);
            ptr_array1GPU = pointer_dict.at(array1GPU);
            ptr_array2Host = pointer_dict.at(array2Host);
            ptr_array2GPU = pointer_dict.at(array2GPU);
        }
        catch(const exception& e) {
            cerr << e.what() << '\n';
            return;
        }

        // Buffer to messages
        char buffer_message [200];
        
        /* Move data from GPU to CPU */
        cudaError_t cuda_err;
        cuda_err = cudaMemcpy(ptr_array1Host, ptr_array1GPU, total_bytes_sun, cudaMemcpyDeviceToHost);
        if (cuda_err != cudaSuccess)
        {
            sprintf(buffer_message, "cudaMemcpy %s <- %s returned error code %d, line(%d)\n", array1Host.c_str(), array1GPU.c_str(), cuda_err, __LINE__);
            logger->print_message(__func__, __LINE__, buffer_message, error_type);
        }

        cuda_err = cudaMemcpy(ptr_array2Host, ptr_array2GPU, total_bytes_py, cudaMemcpyDeviceToHost);
        if (cuda_err != cudaSuccess)
        {
            sprintf(buffer_message, "cudaMemcpy %s <- %s returned error code %d, line(%d)\n", array2Host.c_str(), array2GPU.c_str(), cuda_err, __LINE__);
            logger->print_message(__func__, __LINE__, buffer_message, error_type);
        }

        /* Sum each array and compare the size */

        // Sundials
        // ySun = {T0, Y00, Y01, ... Y0(NSP-1), T1, Y10, Y11, ... Y1(NSP-1), ...}
        double total_array1{0.0f};
        for (int i = 0; i < systems * nsp; i++) {
            total_array1 += ptr_array1Host[i];
        }
        sprintf(buffer_message, "Total sum of %s array = %E", array1Host.c_str(), total_array1);
        logger->print_message(__func__, __LINE__, buffer_message, info_type);

        // Pyjac
        // yPy  = {T0,  T1,  T2, ...  T(nSystems), Y00, Y10, Y20, ..., Y(nSystems)0, Y01, Y11, Y21, ..., Y(nSystems)1, ...}
        double total_array2{0.0f};
        for (int i = 0; i < nsp; i++) {
            for (int j = 0; j < systems; j++) {
                total_array2 += ptr_array2Host[i * padded + j];
            }
        }
        sprintf(buffer_message, "Total sum of %s array = %E", array2Host.c_str(), total_array2);
        logger->print_message(__func__, __LINE__, buffer_message, info_type);

        double error = abs(total_array1 - total_array2);

        sprintf(buffer_message, "Difference between %s and %s = %E", array1Host.c_str(), array2Host.c_str(), error);
        logger->print_message(__func__, __LINE__, buffer_message, info_type);

        if (error > 0.0) {
            sprintf(buffer_message, "The sum of each array (%s and %s) is different, this is the result of subtraction: %E", array1Host.c_str(), array2Host.c_str(), error);
            logger->print_message(__func__, __LINE__, buffer_message, error_type);
        }
    }

    /**********************************************************************/

    /************************* JacSunJacPyjac ***************************/
    /* PUBLIC */

    // Constructor
    JacSunJacPyjac::JacSunJacPyjac(double *_JacSunHost, double *_JacSunGPU, double *_JacPyjacHost, double *_JacPyjacGPU) \
        : JacSunHost(_JacSunHost), JacSunGPU(_JacSunGPU), JacPyjacHost(_JacPyjacHost), JacPyjacGPU(_JacPyjacGPU) {

        logger = make_unique<Utils::Logger>("Jac matrix");
    }


    void JacSunJacPyjac::set_simulation_size(int num_systems, int num_species, int _padded) {
        systems = num_systems;
        nsp = num_species;

        // Padded is the number of systems used by pyjac
        padded = _padded;

        total_bytes_sun = systems * nsp * nsp * sizeof(double);
        total_bytes_py = padded * nsp * nsp * sizeof(double);
    }

    /* Level can be "error" or "info" */
    void JacSunJacPyjac::set_logger_level(string level) {
        logger->set_logger_level(level);
    }

    void JacSunJacPyjac::compare_matrices() {
        // Buffer to messages
        char buffer_message [200];
        
        /* Move data from GPU to CPU */
        cudaError_t cuda_err;
        cuda_err = cudaMemcpy(JacSunHost, JacSunGPU, total_bytes_sun, cudaMemcpyDeviceToHost);
        if (cuda_err != cudaSuccess)
        {
            sprintf(buffer_message, "cudaMemcpy JacSunHost <- JacSunGPU returned error code %d, line(%d)\n", cuda_err, __LINE__);
            logger->print_message(__func__, __LINE__, buffer_message, error_type);
        }

        cuda_err = cudaMemcpy(JacPyjacHost, JacPyjacGPU, total_bytes_py, cudaMemcpyDeviceToHost);
        if (cuda_err != cudaSuccess)
        {
            sprintf(buffer_message, "cudaMemcpy JacPyjacHost <- JacPyjacGPU returned error code %d, line(%d)\n", cuda_err, __LINE__);
            logger->print_message(__func__, __LINE__, buffer_message, error_type);
        }

        /* Sum each matrix and compare the size */
        /* Compare by chunks (individual matrices of 53x53 size) */
        // Sundials
        vector<double> individual_matrix_sum1(systems, 0.0f);
        double total_matrix1(0.0f);
        for (int i = 0; i < systems; i++) {
            for (int j = 0; j < nsp * nsp; j++) {
                individual_matrix_sum1[i] += JacSunHost[i * nsp * nsp + j];
            }
            total_matrix1 += individual_matrix_sum1[i];
        }

        // Pyjac
        // yPy  = {T0,  T1,  T2, ...  T(nSystems), Y00, Y10, Y20, ..., Y(nSystems)0, Y01, Y11, Y21, ..., Y(nSystems)1, ...}
        vector<double> individual_matrix_sum2(systems, 0.0f);
        double total_matrix2(0.0f);
        for (int j = 0; j < systems ; j++) {
            for (int i = 0; i < nsp*nsp; i++) {
                individual_matrix_sum2[j] += JacPyjacHost[i * padded + j];
                
            }
            total_matrix2 += individual_matrix_sum2[j];
        }

        for (int i = 0; i < systems; i++) {
            sprintf(buffer_message, "Sum of matrix number %d of JacSunHost = %E", i, individual_matrix_sum1[i]);
            logger->print_message(__func__, __LINE__, buffer_message, info_type);

            sprintf(buffer_message, "Sum of matrix number %d of JacPyjacHost = %E", i, individual_matrix_sum2[i]);
            logger->print_message(__func__, __LINE__, buffer_message, info_type);

            double individual_error = abs(individual_matrix_sum1[i] - individual_matrix_sum2[i]);
            sprintf(buffer_message, "Difference between JacSunHost and JacPyjacHost matrix number %d = %E", i, individual_error);
            logger->print_message(__func__, __LINE__, buffer_message, info_type);

            if (individual_error > 0.0) {
                sprintf(buffer_message, "The sum of the jacobian matrix number %d (JacSunHost and JacPyjacHost) is different, this is the result of subtraction: %E", i, individual_error);
                logger->print_message(__func__, __LINE__, buffer_message, error_type);
            }
        }


        sprintf(buffer_message, "Total sum of JacSunHost matrix = %E", total_matrix1);
        logger->print_message(__func__, __LINE__, buffer_message, info_type);

        sprintf(buffer_message, "Total sum of JacPyjacHost matrix = %E", total_matrix2);
        logger->print_message(__func__, __LINE__, buffer_message, info_type);

        double error = abs(total_matrix1 - total_matrix2);

        sprintf(buffer_message, "Difference between JacSunHost and JacPyjacHost = %E", error);
        logger->print_message(__func__, __LINE__, buffer_message, info_type);

        if (error > 0.0) {
            sprintf(buffer_message, "The sum of each matrix (JacSunHost and JacPyjacHost) is different, this is the result of subtraction: %E", error);
            logger->print_message(__func__, __LINE__, buffer_message, error_type);
        }
    }
} // namespace Testing

