//#include "ODEIntegrator/Integrators/Drivers/CVodeGPUDriver.cuh"

#include "Mechanism/GPU/header.cuh"

#include "cuda_runtime.h"

#include <iostream>

void init_memory_gpu(int num_systems, mechanism_memory *pyjac_mem) {
    cudaError_t cuda_err;

    /* Allocate Memory */

    cuda_err = cudaMalloc((void **) &(pyjac_mem->y), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->dy), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->conc), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->fwd_rates), FWD_RATES * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->rev_rates), REV_RATES * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->spec_rates), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->cp), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->h), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->dBdT), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->jac), NSP * NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->var), 1 * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->J_nplusjplus), NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMalloc((void **) &(pyjac_mem->pres_mod), PRES_MOD_RATES * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMalloc error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    /* Init memory */

    cuda_err = cudaMemset(pyjac_mem->spec_rates, 0, NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMemset error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMemset(pyjac_mem->dy, 0, NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMemset error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaMemset(pyjac_mem->jac, 0, NSP * NSP * num_systems * sizeof(double));
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaMemset error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }
}

void free_memory_gpu(mechanism_memory *pyjac_mem) {
    cudaError_t cuda_err;

    cuda_err = cudaFree(pyjac_mem->y);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->dy);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->conc);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->fwd_rates);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->rev_rates);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->spec_rates);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->cp);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->h);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->dBdT);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->jac);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->var);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->J_nplusjplus);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

    cuda_err = cudaFree(pyjac_mem->pres_mod);
    if (cuda_err != cudaSuccess) {
        std::cout << "Error: cudaFree error code " << cuda_err << " on function \"";
        std::cout << __func__ << "\" (line: " << __LINE__ << ")" << std::endl;
    }

}



