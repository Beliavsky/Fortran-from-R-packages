#include "header.h"

/* These function estimate the extremal coefficient by using: (a) the
   method proposed by Smith in his unpublished manuscript and (b) the
   schlather and tawn estimator.
*/

void extCoeffSmith(double *frech, int *nObs, int *nSite,
		   double *extCoeff){
  const int nPair = *nSite * (*nSite - 1) / 2;

  for (int currentPair=0;currentPair<nPair;currentPair++){
    int i, j;
    getSiteIndex(currentPair, *nSite, &i, &j);
    for (int k=0;k<*nObs;k++)
      extCoeff[currentPair] +=  fmin2(frech[i * *nObs + k], frech[j * *nObs + k]);

    extCoeff[currentPair] = *nObs / extCoeff[currentPair];
  }

  return;
}

void extCoeffST(double *frech, double *xbar, double *z, double *theta,
		int *nObs, double *dns){

  for (int i=0;i<*nObs;i++){
    double frechScaled = fmax2(frech[i] * xbar[0], frech[i + *nObs] * xbar[1]);

    if (frechScaled > *z)
      *dns +=  -log(*theta) + *theta / frechScaled;

    else
      *dns += *theta / *z;
  }

  return;
}







