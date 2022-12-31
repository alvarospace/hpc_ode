#pragma once

#include <vector>
#include <memory>
#include <string>
#include <sstream>

#include "ODEIntegrator/Integrators/Integrator.hpp"
#include "ODEIntegrator/Integrators/MechanismDriver.hpp"
#include "ODEIntegrator/Mesh/Mesh.hpp"
#include "ODEIntegrator/Integrators/CVodeCPUDataModels.hpp"

class CVodeIntegrator : public Integrator, public MechanismDriver {
    public:
        virtual void init(std::shared_ptr<Context> ctx, IntegratorConfig config) override;
        virtual void integrate(double t0, double t) override;
        virtual dydt_driver dydt_func() override;
        virtual jacobian_driver jacobian_func() override;

        virtual void check_return_value(void* returnValue, std::string const funcName, int const opt);

    protected:
        virtual std::vector<std::vector<double>> data_transfer_from_mesh();
        virtual void data_transfer_to_mesh(std::vector<std::vector<double>> systemsData);
        void integrateSystem(double* system, double dt);
        double last_specie_calculation(double* species);

    private:
        std::unique_ptr<UserData> uData;
};