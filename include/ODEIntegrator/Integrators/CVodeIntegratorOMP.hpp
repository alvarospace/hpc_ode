#pragma once

#include <vector>

#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"

class CVodeIntegratorOMP : public CVodeIntegrator {
    public:
        void integrate(double t0, double t) override;
    
    private:
        std::vector<std::vector<double>> data_transfer_from_mesh() override;
        void data_transfer_to_mesh(std::vector<std::vector<double>> systemsData) override;
};