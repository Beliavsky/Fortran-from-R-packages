#ifndef logsumexp2_h
#define logsumexp2_h

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


void logsumexp2_c(double *a, double *b, int n, double *res);
SEXP logsumexp2(SEXP a, SEXP b, SEXP n);

#endif /* logsumexp2_h */
