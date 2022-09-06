#pragma once

#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Integrators/MechanismDriver.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"

#include <vector>
#include <memory>
#include <string>
#include <sstream>

using std::vector;
using std::string;

class CVodeIntegrator : public Integrator, public MechanismDriver {
    public:
        virtual void init(IntegratorConfig config) override;
        virtual void integrate(double t0, double t) override;
        virtual void clean() override {}
        virtual dydt_driver dydt_func() override;
        virtual jacobian_driver jacobian_func() override;

        virtual void check_return_value(void* returnValue, string const funcName, int const opt);

        struct UserData {
            double pressure;
        };
        

    protected:
        virtual vector<vector<double>> data_transfer_from_mesh(Mesh& mesh);
        virtual void data_transfer_to_mesh(Mesh& mesh, vector<vector<double>> systemsData);
        void integrateSystem(double* system, double dt);
        double last_specie_calculation(double* species);

        std::unique_ptr<UserData> uData;
        
};