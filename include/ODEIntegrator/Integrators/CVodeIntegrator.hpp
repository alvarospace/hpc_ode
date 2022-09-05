#pragma once

#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Integrators/MechanismDriver.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "sundials/sundials_linearsolver.h"

#include <vector>
#include <memory>

using std::vector;

class CVodeIntegrator : public Integrator, public MechanismDriver {
    public:
        virtual void init(IntegratorConfig config) override;
        virtual void integrate(double t0, double t) override;
        virtual void clean() override {}
        virtual dydt_driver dydt_func() override;
        virtual jacobian_driver jacobian_func() override;

        struct UserData {
            double pressure;
        };
        

    private:
        vector<vector<double>> data_transfer_from_mesh(Mesh& mesh);
        void data_transfer_to_mesh(Mesh& mesh, vector<vector<double>> systemsData);
        void integrateSystem(double* system, double dt);
        double last_specie_calculation(double* species);

        std::unique_ptr<UserData> uData;
        
};