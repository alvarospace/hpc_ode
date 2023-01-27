#ifndef MECHANISM_cuh
#define MECHANISM_cuh

#ifdef __GNUG__
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "launch_bounds.cuh"
#include "gpu_macros.cuh"
#endif

struct mechanism_memory {
  double * y;
  double * dy;
  double * conc;
  double * fwd_rates;
  double * rev_rates;
  double * spec_rates;
  double * cp;
  double * h;
  double * dBdT;
  double * jac;
  double * var;
  double * J_nplusjplus;
  double * pres_mod;
};

//last_spec 52
/* Species Indexes
0  H2
1  H
2  O
3  O2
4  OH
5  H2O
6  HO2
7  H2O2
8  C
9  CH
10  CH2
11  CH2(S)
12  CH3
13  CH4
14  CO
15  CO2
16  HCO
17  CH2O
18  CH2OH
19  CH3O
20  CH3OH
21  C2H
22  C2H2
23  C2H3
24  C2H4
25  C2H5
26  C2H6
27  HCCO
28  CH2CO
29  HCCOH
30  N
31  NH
32  NH2
33  NH3
34  NNH
35  NO
36  NO2
37  N2O
38  HNO
39  CN
40  HCN
41  H2CN
42  HCNN
43  HCNO
44  HOCN
45  HNCO
46  NCO
47  N2
48  AR
49  C3H7
50  C3H8
51  CH2CHO
52  CH3CHO
*/

//Number of species
#define NSP 53
//Number of variables. NN = NSP + 1 (temperature)
#define NN 54
//Number of forward reactions
#define FWD_RATES 325
//Number of reversible reactions
#define REV_RATES 309
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 50

//Must be implemented by user on a per mechanism basis in mechanism.cu
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

