// This program solves equations of type M y'= f(y) 

#include "cantera/zerodim.h"
#include "cantera/IdealGasMix.h"

void run()
{
    // Size of the problem Y (Unknowns)
    int nsp;
    // Input Constants 
    double p = 101325;
    double temp = 1000;
    // Final output 
    double tempf
    // Final time
    double tf   = 1e-3;

    // Objects (need to be made into pointers and vectorised)
    sol = Cantera::newSolution("gri30.xml", "gri30", "None")
    Cantera::IdealGasConstPressureReactor r;
    Cantera::ReactorNet net;

    // Get size of unknowns 
    nsp = sol.nSpecies();
    // Create Vector of unknowns
    Cantera::vector_fp x(nsp, 0.0);

    // Define initial Values of unknows (Corresponds to phi = 1.0, perfect fuel air mixture) 
    x[gas[0]->speciesIndex("CH4")] = 0.19011;
    x[gas[0]->speciesIndex("O2")] = 0.095057;
    x[gas[0]->speciesIndex("N2")] = 0.71483;

    // Set initial Values of unknows in vector 
    for (int i=0; i<n_size; i++) {
        sol.setState_TPY(temp,p,x.data());
    }


    // operate on pointers
    r.insert(sol);
    net.addReactor(r);

    // Advance Solution by time tf
    net.advance(tf);

    // Solution
    tempf =  r.temperature();
    std::cout << "Final Temperature "<< tempf << std::endl;

}

int main()
{
    try {
        run();
        return 0;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        std::cout << " terminating... " << std::endl;
        appdelete();
        return 1;
    }
}
