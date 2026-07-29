#ifndef eval_gdenxmis_integrand_h
#define eval_gdenxmis_integrand_h

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


void eval_gdensmix_integrand_c(double *x, double *mix, int *groupings,
                               double *factorinv, int d, int N, int n,
                               double lconst, int returnall, double *res);

SEXP eval_gdensmix_integrand_returnall(SEXP x, SEXP mix, SEXP groupings, SEXP factorinv,
                                       SEXP d, SEXP N, SEXP n, SEXP lconst);

SEXP eval_gdensmix_integrand_LSE(SEXP x, SEXP mix, SEXP groupings, SEXP factorinv,
                                 SEXP d, SEXP N, SEXP n, SEXP lconst);


#endif /* eval_gdenxmis_integrand_h */
