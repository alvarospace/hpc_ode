# CVODES GPU HINTS 

Initializing y (unknown) =   

int neq_tot = (NUM_SPECIES + 1) * ncells;
y = N_VNew_Serial(neq_tot, The_Sundials_Context());


Options for linear solver 

1. Cuda based liner solver from sundials = SUNLinSol_cuSolverSp_batchQR

2. Magma solver from SUNLinSol_MagmaDense






