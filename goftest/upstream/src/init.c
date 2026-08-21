/* 
   Native symbol registration table for goftest package

*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

void ADprobExactInf(double *, int *, double *);
void ADprobN(double *, int *, int *, double *);
void ADprobApproxInf(double *, int *, double *);
void ADtestR(double *, int *, double *, double *);

static const R_CMethodDef CEntries[] = {
    {"ADprobExactInf",   (DL_FUNC) &ADprobExactInf,   3},
    {"ADprobN",          (DL_FUNC) &ADprobN,          4},
    {"ADprobApproxInf",  (DL_FUNC) &ADprobApproxInf,  3},
    {"ADtestR",          (DL_FUNC) &ADtestR,          4},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {NULL, NULL, 0}
};

void R_init_goftest(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
  
