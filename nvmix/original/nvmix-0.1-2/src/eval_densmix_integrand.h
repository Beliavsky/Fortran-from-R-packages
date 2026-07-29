#ifndef eval_densmix_integrand_h
#define eval_densmix_integrand_h

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


void eval_densmix_integrand_c(double *W, double *maha2_2, int current_n, int n,
                             int d, int k, double *lconst, double *ldensities,
                             double *c);
SEXP eval_densmix_integrand(SEXP W, SEXP maha2_2, SEXP current_n, SEXP n, SEXP d,
                           SEXP k, SEXP lconst);

#endif /* eval_densmix_integrand_h */
