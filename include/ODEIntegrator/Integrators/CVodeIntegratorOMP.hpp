#pragma once

#include "ODEIntegrator/Integrators/CVodeIntegrator.hpp"

class CVodeIntegratorOMP : public CVodeIntegrator {
    public:
        void integrate(double t0, double t) override;
    
    private:
        vector<vector<double>> data_transfer_from_mesh(Mesh& mesh) override;
        void data_transfer_to_mesh(Mesh& mesh, vector<vector<double>> systemsData) override;
};