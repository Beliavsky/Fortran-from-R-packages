#ifndef R_GMM_H
#define R_GMM_H

#include <R_ext/RS.h>


void F77_SUB(wu)(double *gt, double *tol, int *maxit,
		  int *n, int *q, double *k, int *conv,
		  double *obj, double *lam);
		 
void F77_SUB(lamcuep)(double *gt, int *n, int *q, double *k,
		       int *maxit, int *conv, double *lam,
		       double *pt, double *obj);

#endif		     
