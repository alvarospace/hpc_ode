# HPC_Cvode

Implementation of stand-alone CVODE integrators for CPU and GPU implementations

This implementation of the integrator needs functions generated from an external library called PyJac. Pyjac containts information about the chemisty, this is intended to replace cantera for chemistry data. 

## Dependencies 

main.cpp - Contains csv reading funciton and integrators from sundials (6.1) 

cvodes_dydt.c - Wrapper class for cvode to get dy/dt
 
cvodes_jac.c - Wrapper class for cvode to get Jacobian 

writeLibrary.py - A python script to compile sources generated from pyjac, complie the out/ directory 

out  - directory that contains source and header files generated from pyjac for gri30 chemical mechanism


