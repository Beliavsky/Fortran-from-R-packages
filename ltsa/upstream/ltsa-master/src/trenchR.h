#ifndef _TRENCH_H_
#define _TRENCH_H_

#include "nrutil.h"
#include <string.h>

#include <R.h>

int durlev (double *c,int n,double *fi,double *v,double *vk,double EPS);
int trenchInv(double *r,int n,double **b,double *v,double EPS);
void upperWedge(double *phi,double sigsq,int n,double **b);
void fromWedgeStorage(int n,double **B);
double trenchDet(double *r,int n,double *v);
double dot(int n,double* u,double* v);
void vecmat(int nc,double *v, double **m,double *u);

#endif /* _TRENCH_H_ */
