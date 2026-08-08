#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "n1qn1.h"

extern void n1qn1_(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
		   int mode[], int niter[], int nsim[], int imp[], double zm[], int izs[],
		   float rzs[], double dzs[], int id[]);
void n1qn1F(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
		   int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[],
		   float rzs[], double dzs[]) {
  int id = 0;
  n1qn1_(simul, n, x, f, g, var, eps, mode, niter, nsim, imp, zm, izs, rzs, dzs, &id);
}

void n1qn1F2(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
		   int mode[], int niter[], int nsim[], int imp[], double zm[], int izs[],
		   float rzs[], double dzs[]) {
  int id = 0;
  n1qn1_(simul, n, x, f, g, var, eps, mode, niter, nsim, imp, zm, izs, rzs, dzs, &id);
}
/* .C calls */
extern SEXP n1qn1_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _n1qn1_ptr(void) {
  int pro = 0;  // Counter for the number of PROTECT calls

  // Create an external pointer
  SEXP n1qn1c_n1qn1F = PROTECT(R_MakeExternalPtrFn((DL_FUNC)&n1qn1F, R_NilValue, R_NilValue)); pro++;

  SEXP n1qn1c_n1qn1F2 = PROTECT(R_MakeExternalPtrFn((DL_FUNC)&n1qn1F2, R_NilValue, R_NilValue)); pro++;

  SEXP n1qn1c_n1qn1_ = PROTECT(R_MakeExternalPtrFn((DL_FUNC)&n1qn1_, R_NilValue, R_NilValue)); pro++;

#define nVec 3

  SEXP ret = PROTECT(Rf_allocVector(VECSXP, nVec)); pro++;
  SEXP retN = PROTECT(Rf_allocVector(STRSXP, nVec)); pro++;

  SET_VECTOR_ELT(ret, 0, n1qn1c_n1qn1F);
  SET_STRING_ELT(retN, 0, Rf_mkChar("n1qn1F"));

  SET_VECTOR_ELT(ret, 1, n1qn1c_n1qn1F2);
  SET_STRING_ELT(retN, 1, Rf_mkChar("n1qn1F2"));

  SET_VECTOR_ELT(ret, 2, n1qn1c_n1qn1_);
  SET_STRING_ELT(retN, 2, Rf_mkChar("n1qn1_"));

#undef nVec

  // Set the names attribute of the list
  Rf_setAttrib(ret, R_NamesSymbol, retN);

  // Unprotect all protected objects
  UNPROTECT(pro);

  // Return the list of external pointers
  return ret;

}

void R_init_n1qn1(DllInfo *dll)
{
  R_CallMethodDef callMethods[]  = {
    {"_n1qn1_ptr", (DL_FUNC) &_n1qn1_ptr, 0},
    {"n1qn1_wrap", (DL_FUNC) &n1qn1_wrap, 13},
    {NULL, NULL, 0}
  };
  R_RegisterCCallable("n1qn1","n1qn1F", (DL_FUNC) &n1qn1F);
  R_RegisterCCallable("n1qn1","n1qn1F2", (DL_FUNC) &n1qn1F2);
  R_RegisterCCallable("n1qn1","n1qn1_", (DL_FUNC) &n1qn1_);
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
