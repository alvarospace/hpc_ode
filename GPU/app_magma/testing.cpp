#include <iostream>
#include <string>
#include <app_magma/testing.hpp>

using namespace std;

namespace Testing
{
    
    void mesh_vs_ysun(const std::shared_ptr<Utils::ThermData>& mesh, const double *ysun) {

        Utils::Logger infoLogger("Testing", "INFO");
        double error{0.0f};
        
        string message = "Accumulated error between data = " + to_string(error);
        infoLogger.print_message(__func__, message);

    }
} // namespace Testing

