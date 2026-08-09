#include<R_ext/Rdynload.h>
#ifndef R_R_H
#  include <R.h>
#endif


void F77_NAME(ldp)(double*, double*, int*, int*, int*,  double*, 
      double*, double*, int*, int*, int*, int*, int*);
     
void F77_NAME(lsei)(int*, int*, int*, int*, double*, double*,
      double*, double*, double*, double*, double*, int*, int*, 
      int*, int*, double*, double*, int*, double*, int*, int*);

void F77_NAME(xnnls)(double*, int*, int*, int*, double*, 
      double*, double*, double*, double*, int*, int*, int*);             

void F77_NAME(dgbsv)(int*, int*, int*, int*, double*, int*,
      int*, double*, int*, int*);

void F77_NAME(block)(int*, double*, int*, int*, double*, int*, 
      int*, int*, double*, int*, int*, int*,
      int*, double*, double*, int*, double*, double*);
	  
void F77_NAME(dgtsv)(int*, int*, double*,
      double*, double*, double*, int*, int*);
          
R_FortranMethodDef fortranMethods[] = {
 {"ldp",(DL_FUNC) &F77_SUB(ldp), 13},
 {"lsei", (DL_FUNC) &F77_SUB(lsei), 21},
 {"xnnls", (DL_FUNC) &F77_SUB(xnnls), 12},
 {"dgbsv", (DL_FUNC) &F77_SUB(dgbsv), 10},
 {"block", (DL_FUNC) &F77_SUB(block), 18},
 {"dgtsv", (DL_FUNC) &F77_SUB(dgtsv), 8},
 {NULL, NULL, 0}
};
void R_init_limSolve(DllInfo *info) {
  R_registerRoutines(info, NULL, NULL, fortranMethods, NULL);
  R_useDynamicSymbols(info, FALSE); // disable dynamic searching  
}
