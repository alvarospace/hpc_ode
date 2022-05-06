#include <cvodes/cvodes.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>

extern "C"{
#include "dydt.h"
#include "jacob.h"
}

#include <ctime>
#include <iostream>
#include <stdio.h>
#include <fstream>

#include <string>
#include <vector>
#include <sstream>

using std::ofstream;
using namespace std;


int eval_jacob_cvodes(realtype t, N_Vector y, N_Vector ydot, SUNMatrix Jac, void* f, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int dydt_cvodes(realtype t, N_Vector y, N_Vector ydot, void* f);

int main()
{
   ofstream file2;	
   // Sundials declarations  
   void *cvode_mem = NULL;
   SUNContext sunctx;
   SUNMatrix A;
   SUNLinearSolver LS;
   N_Vector y, ydot;
   realtype t;
   // Tolorences 
   realtype reltol = RCONST(1.0e-6);
   realtype abstol = RCONST(1.0e-10);
   int m_maxsteps = 10000;
   
   A = NULL;
   LS = NULL;
   y = NULL;
   double Ytest[NSP] = {0};
   ydot = NULL;
   y = N_VNew_Serial(NSP,sunctx);
   ydot = N_VNew_Serial(NSP,sunctx);
   int retval = SUNContext_Create(NULL, &sunctx);
   
   double tfinal;
   double f;
   int cante;
   cante = 0;
   //cante = 1;

   // Varibales for reading file 
   vector<vector<double>> content;
   vector<double> row;

   cout << " tfinal " << endl; 
   cin >> tfinal;
   
   // Read File
   string fname="../res-est.csv";
   string line, word;
   double data;
   fstream file (fname, ios::in);
   while(getline(file, line)) {
        row.clear();
        stringstream str(line);
        while(getline(str, word, ',')){
             data = stod(word);
             row.push_back(data);
	}
        content.push_back(row);
   }
   // Pressure fixed to 1 bar
   f = 101325.15;
   string fnameout = "./time_" + std::to_string(tfinal) + ".txt";
   file2.open(fnameout,std::ios_base::app);

   for(int i=0;i<content.size();i++){	   
      // Pyjac + CVODE integration     
      NV_Ith_S(y,0) = content[i][NSP];
      for (int j=0 ; j < NSP-1 ; j++){
          NV_Ith_S(y,j+1) = content[i][j];    
      }
      cvode_mem = CVodeCreate(CV_BDF,sunctx);
      int flag = CVodeInit(cvode_mem, dydt_cvodes, 0, y);
      flag = CVodeSStolerances(cvode_mem, reltol, abstol);
      A = SUNDenseMatrix(NSP, NSP,sunctx);
      LS = SUNLinSol_Dense(y, A,sunctx);
      flag = CVodeSetLinearSolver(cvode_mem, LS, A);
      flag = CVodeSetJacFn(cvode_mem,  eval_jacob_cvodes);
      flag = CVodeSetMaxNumSteps(cvode_mem, m_maxsteps);
      flag = CVodeSetUserData(cvode_mem, &f);
      t = 0.0;
      clock_t c_start = clock();
      flag = CVode(cvode_mem, tfinal, y, &t, CV_NORMAL);
      clock_t c_end = clock();
      double timee = c_end - c_start;
      file2 << timee << std::endl;
   }
   file2.close();
   N_VDestroy(y);
   CVodeFree(&cvode_mem);
   return 0;
}

int eval_jacob_cvodes(double t, N_Vector y, N_Vector ydot, SUNMatrix jac, void* f, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{       
	double* local_y = NV_DATA_S(y);
	double temp[NSP][NSP] = {0};
	eval_jacob((double)t, *(double*)f, local_y, (double*)temp);

	for (int i=0; i < NSP ; i++) {
	    for (int j=0; j < NSP ; j++) {
	        SM_ELEMENT_D(jac,j,i) = temp[i][j];
	    }
	}

	return 0;
}

int dydt_cvodes(realtype t, N_Vector y, N_Vector ydot, void* f)
{
	double* local_y = NV_DATA_S(y);
	double* local_dy = NV_DATA_S(ydot);
	dydt((double)t, *(double*)f, local_y, local_dy);
	return 0;
}
