#ifndef precond_h
#define precond_h

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>



void precond(double *lower, double *upper, double *scale,
             double *C, double *meansqrtmix, double *tol,
             int *d_, int *perm, int *status);
void swap(double *lower, double *upper, double *meansqrtmix, int *perm,
          double *scale, int d, int i, int j);

#endif /* precond_h */




