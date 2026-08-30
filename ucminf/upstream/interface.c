/*
  Stig Bousgaard Mortensen
  DTU Informatics
  sbm@imm.dtu.dk

  Based on R-package 'FortranCallsR' by
  Diethelm Wuertz, ETH Zurich, www.rmetrics.org

  Modifications by Douglas Bates <bates@stat.wisc.edu>, Nov. 2010
  Modifications by Tomas Kalibera <tomas.kalibera@gmail.com> Aug. 2016.
  Modifications by K Hervé Dakpo <k-herve.dakpo@inrae.fr> Dec. 2023.
*/


#include <R.h>
#include <Rinternals.h>  //R internal structures
#include <R_ext/RS.h>    //F77_CALL etc.
#include <Rversion.h>

// Declare FORTRAN routine for use in C
extern void F77_NAME(ucminf)(int*, double[], double*, double[],
			     int*,double[],int*,int*,int*,double[],SEXP);

/*-------------------------------------------------------------------------------
  Define C functions that calls user defined function in R
*/

/*
 Compatibility helper:
 Rf_findVarInFrame/findVarInFrame are non-API entry points.
 Use R_getVar(..., inherits = FALSE) on R >= 4.5.0 and fall back on
 older R releases.

 Values looked up in 'rho' do not need PROTECT here: 'rho' is an argument
 to .Call and therefore remains reachable/protected for the duration of the
 native call.
 */
static SEXP getVarInFrame(SEXP rho, const char *name)
{
  SEXP sym = install(name);
#if R_VERSION >= R_Version(4, 5, 0)
  return R_getVar(sym, rho, FALSE);
#else
  return findVarInFrame(rho, sym);
#endif
}

void installPar(int nn, double x[], SEXP rho) {
    int i;
    SEXP PAR = getVarInFrame(rho, ".x");
    double *xpt = REAL(PAR);
    if (LENGTH(PAR) != nn)
	error("Dimension mismatch, length(.x) = %d != n = %d", LENGTH(PAR), nn);
    for (i = 0; i < nn; i++) xpt[i] = x[i] ;
}

void F77_SUB(func)(int *n, double x[], double *value, SEXP rho) {
  SEXP dotf, ans;

  installPar(*n, x, rho);
  dotf = getVarInFrame(rho, ".f");
  PROTECT(ans = eval(dotf, rho));
  *value = asReal(ans);
  UNPROTECT(1);
}

void F77_SUB(usrgr)(int *n, double x[], double grval[], SEXP rho) {
  SEXP dotgr, OUT;
  int i, nn = *n;
  double *grv;

  installPar(nn, x, rho);
  dotgr = getVarInFrame(rho, ".gr");
  PROTECT(OUT = eval(dotgr, rho));
  if (LENGTH(OUT) != nn || !isReal(OUT))
    error("gradient evaluation must return a numeric vector of length %d", nn);
  grv = REAL(OUT);
  for (i = 0; i < nn; i++) grval[i] = grv[i];
  UNPROTECT(1);
}

SEXP mfopt(SEXP rho) {
  int n    = asInteger(getVarInFrame(rho, ".n"));
  int iw   = asInteger(getVarInFrame(rho, ".iw"));
  int grad = asInteger(getVarInFrame(rho, ".grad"));

  SEXP EPS    = getVarInFrame(rho, ".eps");
  SEXP GRSTEP = getVarInFrame(rho, ".grstep");
  SEXP PAR    = getVarInFrame(rho, ".par");
  SEXP icontr = getVarInFrame(rho, ".icontr");
  SEXP maxfun = getVarInFrame(rho, ".maxfun");
  SEXP dx     = getVarInFrame(rho, ".stepmax");
  SEXP W      = getVarInFrame(rho, ".w");

  if (LENGTH(EPS) < 2 || !isReal(EPS))
    error(".eps must be a numeric vector of length >= 2");
  if (LENGTH(GRSTEP) < 2 || !isReal(GRSTEP))
    error(".eps must be a numeric vector of length >= 2");
  if (LENGTH(PAR) != n || !isReal(PAR))
    error("Dimension mismatch, length(.par) = %d != n = %d", LENGTH(PAR), n);
  if (LENGTH(W) != iw || !isReal(W))
    error("Dimension mismatch, length(.w) = %d != .iw = %d", LENGTH(W), iw);

  /* duplicate dx, maxfun, .w because they are input/output arguments */
  PROTECT(maxfun = duplicate(maxfun));
  defineVar(install(".maxfun"), maxfun, rho);
  PROTECT(dx = duplicate(dx));
  defineVar(install(".stepmax"), dx, rho);
  PROTECT(W = duplicate(W));
  defineVar(install(".w"), W, rho);
  UNPROTECT(3); /* now protected via rho */

  F77_CALL(ucminf)(&n, REAL(PAR), REAL(dx), REAL(EPS), INTEGER(maxfun), REAL(W),
           &iw, INTEGER(icontr), &grad, REAL(GRSTEP), rho);
  return R_NilValue;
}

void F77_SUB(prtrac)(int *neval, double *fx, double *nmg, int *n, double x[]) {
  int i, nn = *n;
  Rprintf(" neval = %3d, F(x) =%11.4e, max|g(x)| =%11.4e\n", *neval, *fx, *nmg);
  Rprintf(" x =%11.4e", x[0]);
  for (i = 1; i < nn; i++) Rprintf(",%11.4e", x[i]);
  Rprintf("\n");
}

void F77_SUB(prline)(double *a, double sl[]) {
  Rprintf(" Line search: alpha =%11.4e, dphi(0) =%11.4e, dphi(1) =%11.4e\n",
          *a, sl[0], sl[1]);
}

void F77_SUB(prconv)(void) {
  Rprintf(" Optimization has converged\n");
}

void F77_SUB(prfail)(int *neval) {
  Rprintf(" Optimization stopped after %d function evaluations\n", *neval);
}
