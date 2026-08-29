#include "header.h"

void madogram(double *data, int *nObs, int *nSite, double *mado){
  /* This function computes the madogram */
  const int nPair = *nSite * (*nSite - 1) / 2;

#pragma omp parallel for
  for (int currentPair=0;currentPair<nPair;currentPair++){
    int i,j;
    getSiteIndex(currentPair, *nSite, &i, &j);

      int availableData = 0;
      for (int k=0;k<*nObs;k++){
	if (R_FINITE(data[i * *nObs + k]) && R_FINITE(data[j * *nObs + k])){
	  mado[currentPair] += fabs(data[i * *nObs + k] - data[j * *nObs + k]);
	  availableData++;
	}
      }

      mado[currentPair] *= 0.5 / (double) availableData;
    }

  return;
}

void variogram(double *data, int *nObs, int *nSite, double *vario){
  /* This function computes the (semi) variogram */
  const int nPair = *nSite * (*nSite - 1) / 2;

#pragma omp parallel for
  for (int currentPair=0;currentPair<nPair;currentPair++){
    int i,j;
    getSiteIndex(currentPair, *nSite, &i, &j);

    int availableData = 0;
    for (int k=0; k<*nObs;k++){
      if (R_FINITE(data[i * *nObs + k]) && R_FINITE(data[j * *nObs + k])){
	double dummy = data[i * *nObs + k] - data[j * *nObs + k];
	vario[currentPair] += dummy * dummy;
	availableData++;
      }
    }

    vario[currentPair] *= 0.5 / (double) availableData;
    }

  return;
}

void lmadogram(double *data, int *nObs, int *nSite, double *lambda,
	       int *nLambda, double *lmado){
  /* This function computes the lambda-madogram */
  const int nPair = *nSite * (*nSite - 1) / 2;
  const double cst = 0.5 / *nObs;

#pragma omp parallel for
  for (int currentPair=0;currentPair<nPair;currentPair++){
    int i,j;
    getSiteIndex(currentPair, *nSite, &i, &j);

    for (int k=0;k<*nLambda;k++){
	for (int l=0;l<*nObs;l++){
	  lmado[currentPair * *nLambda + k] +=
	    fabs(R_pow(data[i * *nObs + l], lambda[k]) -
		 R_pow(data[j * *nObs + l], 1 - lambda[k])) -
	    lambda[k] * (1 - R_pow(data[i * *nObs + l], lambda[k])) -
	    (1 - lambda[k]) * (1 - R_pow(data[j * *nObs + l],
					 1 - lambda[k]));
	}

	lmado[currentPair * *nLambda + k] *= cst;
	lmado[currentPair * *nLambda + k] += 0.5 *
	  (1 - lambda[k] + lambda[k] * lambda[k]) / (2 - lambda[k]) /
	  (1 + lambda[k]);
    }
  }

  return;
}
