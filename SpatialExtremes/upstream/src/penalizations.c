#include "header.h"

double penalization(double *penmat, double *beta, double pencoeff,
		    int n, int nppar){
  //This function computes the penalization quantity from the penalty
  //matrix, the penalty coefficient and the vector of coefficient.

  int i,j;
  double penalty = .0;

  for (i=nppar;i<n;i++)
    for (j=nppar;j<n;j++)
      penalty = penalty + beta[j] * penmat[j + i * n] * beta[i];

  return(pencoeff * penalty);
}

double penalization2(double *penmat, double *beta, double pencoeff,
		     int n, int nppar){
  //This function computes the penalization quantity from the penalty
  //matrix, the penalty coefficient and the vector of coefficient.

  int i;
  double penalty = .0;

  for (i=nppar;i<n;i++)
    penalty = penalty + R_pow_di(beta[i], 2);

  return(pencoeff * penalty);
}
