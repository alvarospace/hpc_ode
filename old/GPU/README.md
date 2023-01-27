# CVODES GPU HINTS 

Initializing y (unknown) =   
```
int neq_tot = (NUM_SPECIES + 1) * ncells;
y = N_VNew_Serial(neq_tot, The_Sundials_Context());
```

Options for linear solver 

1. Cuda based liner solver from sundials = SUNLinSol_cuSolverSp_batchQR

2. Magma solver from SUNLinSol_MagmaDense


## Dependencies 

- CVode library from Sundials (6.2)

- PyJac (CUDA version) library for right-hand function, dy/dt = f(t,y), and Analytical Jacobian matrix evaluation

- MAGMA library for batched GPU matrix calculations

## Directories

- "./app_magma/" GPU implementation. It depends from "utils.hpp/cpp" for reading, writing, mesh structures, etc. And "cvode_user.cuh/cu" to link CVode and PyJac libraries

## Build System

Compilation of executables and PyJac library is carried out by CMake.

### Instructions

    - mkdir build
    - cd build
    - cmake ..
    - make -j 4

### Setting CMake to build with atypical Sundials path

    - cmake -DSUNDIALS_EXTERN_ENABLE=ON -DSUNDIALS_EXTERN_LIBS=<path/sundials/cmake/file> ..







