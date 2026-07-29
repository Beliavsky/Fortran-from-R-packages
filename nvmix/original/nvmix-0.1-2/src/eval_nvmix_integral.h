#ifndef eval_nvmix_integral_h
#define eval_nvmix_integral_h

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


void eval_nvmix_integral_c(double *lower, double *upper, int *grouping, int numgroups,
                           double *U, double *rtW, double *rtWant, int n, int d, int r,
                           int *kfactor, double *factor, double ZERO, double ONE,
                           int doant, double *res);

SEXP eval_nvmix_integral(SEXP lower, SEXP upper, SEXP grouping, SEXP numgroups, SEXP U,
                         SEXP rtW, SEXP rtWant, SEXP n, SEXP d, SEXP r, SEXP kfactor,
                         SEXP factor, SEXP ZERO, SEXP ONE, SEXP doant_);

#endif /* eval_nvmix_integral_h */


