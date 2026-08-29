#include "header.h"

/* These functions estimate the spatial dependence parameters by using
   a least square optimization. More precisely, it compares the
   expected extremal coefficients to the "observed" ones.
*/

void fitcovmat2d(double *cov11, double *cov12, double *cov22,
		 int *nPairs, double *distVec, double *extcoeff,
		 double *weights, double *ans){

  double dummy = 0.0,
    *mahalDist = malloc(*nPairs * sizeof(double));

  *ans = - mahalDistFct(distVec, *nPairs, cov11, cov12, cov22,
			mahalDist);

  if (*ans != 0.0){
    *ans = 1e50;
    return;
  }

#pragma omp parallel for reduction(+:dummy)
  for (int i=0;i<*nPairs;i++){
    double res = 2 * pnorm(0.5 * mahalDist[i], 0, 1, 1, 0) - extcoeff[i];
    dummy += res * res / (weights[i] * weights[i]);
  }

  *ans = dummy;

  free(mahalDist);
  return;
}


void fitcovmat3d(double *cov11, double *cov12, double *cov13,
		 double *cov22, double *cov23, double *cov33,
		 int *nPairs, double *distVec, double *extcoeff,
		 double *weights, double *ans){

  double dummy = 0.0,
    *mahalDist = malloc(*nPairs * sizeof(double));

  *ans = - mahalDistFct3d(distVec, *nPairs, cov11, cov12, cov13,
			  cov22, cov23, cov33, mahalDist);

  if (*ans != 0.0)
    return;

#pragma omp parallel for reduction(+:dummy)
  for (int i=0;i<*nPairs;i++){
    double res = 2 * pnorm(0.5 * mahalDist[i], 0, 1, 1, 0) - extcoeff[i];
    dummy += res * res / (weights[i] * weights[i]);
  }

  *ans = dummy;

  free(mahalDist);
  return;
}

void fitcovariance(int *covmod, double *nugget, double *range, double *smooth,
		   double *smooth2, int *nPairs, int *dim, double *dist,
		   double *extcoeff, double *weights, double *ans){

  if (*nugget >= 1){
    *ans = - *nugget * *nugget * MINF;
    return;
  }

  double *rho = malloc(*nPairs * sizeof(double)),
    dummy = 0.0;

  switch (*covmod){
  case 1:
    *ans = -whittleMatern(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 2:
    *ans = -cauchy(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 3:
    *ans = -powerExp(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 4:
    *ans = -bessel(dist, *nPairs, *dim, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 5:
    *ans = -caugen(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, *smooth2, rho);
    break;
  }

  if (*ans != 0.0){
    free(rho);
    return;
  }

#pragma omp parallel for reduction(+:dummy)
  for (int i=0;i<*nPairs;i++){
    double res = 1 + sqrt(0.5 - 0.5 * rho[i]) - extcoeff[i];
    dummy += res * res / (weights[i] * weights[i]);
  }

  *ans = dummy;

  free(rho);
  return;
}

void fittcovariance(int *covmod, double *nugget, double *range, double *smooth,
		    double *smooth2, double *DoF, int *nPairs, int *dim, double *dist,
		    double *extcoeff, double *weights, double *ans){

   if (*nugget >= 1){
    *ans = - *nugget * *nugget * MINF;
    return;
  }

  if (*DoF <= 0){
    *ans = - (1 - *DoF) * (1 - *DoF) * MINF;
    return;
  }

  double dummy = 0.0,
    *rho = malloc(*nPairs * sizeof(double));

  switch (*covmod){
  case 1:
    *ans = -whittleMatern(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 2:
    *ans = -cauchy(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 3:
    *ans = -powerExp(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 4:
    *ans = -bessel(dist, *nPairs, *dim, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 5:
    *ans = -caugen(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, *smooth2, rho);
    break;
  }

  if (*ans != 0.0){
    free(rho);
    return;
  }

#pragma omp parallel for reduction(+:dummy)
  for (int i=0;i<*nPairs;i++){
    double res = 2 * pt(sqrt((1 - rho[i]) * (*DoF + 1) / (1 + rho[i])), *DoF + 1, 1, 0) - extcoeff[i];
    dummy += res * res / (weights[i] * weights[i]);
  }

  *ans = dummy;

  free(rho);
  return;
}

void fiticovariance(int *covmod, double *alpha, double *nugget, double *range,
		    double *smooth, double *smooth2, int *nPairs, int *dim,
		    double *dist, double *extcoeff, double *weights, double *ans){
  /* This computes the least squares for the independent Schlather model */

  if (*alpha > 1){
    *ans = - *alpha * *alpha * MINF;
    return;
  }

  if (*alpha < 0){
    *ans = - (1 - *alpha) * (1 - *alpha) * MINF;
    return;
  }

  if (*nugget >= 1){
    *ans = - *nugget * *nugget * MINF;
    return;
  }

  double dummy = 0.0,
    *rho = malloc(*nPairs * sizeof(double));

  switch (*covmod){
  case 1:
    *ans = -whittleMatern(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 2:
    *ans = -cauchy(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 3:
    *ans = -powerExp(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 4:
    *ans = -bessel(dist, *nPairs, *dim, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 5:
    *ans = -caugen(dist, *nPairs, *nugget, 1 - *nugget, *range, *smooth, *smooth2, rho);
    break;
  }

  if (*ans != 0.0){
    free(rho);
    return;
  }

#pragma omp parallel for reduction(+:dummy)
  for (int i=0;i<*nPairs;i++){
    double res = 2 * *alpha + (1 - *alpha) * (1 + sqrt(0.5 - 0.5 * rho[i])) - extcoeff[i];
    dummy +=  res * res / (weights[i] * weights[i]);
  }

  *ans = dummy;

  free(rho);
  return;
}

void fitgcovariance(int *covmod, double *sigma2, double *sigma2Bound, double *nugget,
		    double *range, double *smooth, double *smooth2, int *nPairs,
		    int *dim, double *dist, double *extcoeff, double *weights,
		    double *ans){
  /* This computes the least squares for the geometric Gaussian model */

  if (*nugget >= 1){
    *ans = - *nugget * *nugget * MINF;
    return;
  }

  double dummy = 0.0,
    *rho = malloc(*nPairs * sizeof(double));

  *ans = -geomCovariance(dist, *nPairs, *dim, *covmod, *sigma2, *sigma2Bound,
			 *nugget, *range, *smooth, *smooth2, rho);

  if (*ans != 0.0){
    free(rho);
    return;
  }

#pragma omp parallel for reduction(+:dummy)
  for (int i=0;i<*nPairs;i++){
    double res = 2 * pnorm(0.5 * rho[i], 0.0, 1.0, 1, 0) - extcoeff[i];
    dummy += res * res / (weights[i] * weights[i]);
  }

  *ans = dummy;

  free(rho);
  return;
}

void fitbrcovariance(double *range, double *smooth, int *nPairs,
		     double *dist, double *extcoeff, double *weights,
		     double *ans){
  /* This computes the least squares for the Brown-Resnick model */

  double dummy = 0.0,
    *rho = malloc(*nPairs * sizeof(double));

  *ans = -brownResnick(dist, *nPairs, *range, *smooth, rho);

  if (*ans != 0.0){
    free(rho);
    return;
  }

  #pragma omp parallel for reduction(+:dummy)
  for (int i=0;i<*nPairs;i++){
    double res = 2 * pnorm(0.5 * rho[i], 0.0, 1.0, 1, 0) - extcoeff[i];
    dummy += res * res / (weights[i] * weights[i]);
  }

  *ans = dummy;
  free(rho);
  return;
}
