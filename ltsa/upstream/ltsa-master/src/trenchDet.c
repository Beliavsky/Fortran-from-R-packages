// Determinant computation of a Toeplitz matrix

// Takes arguments r,n and v 
// Returns logDet, the logarithm of determinant

#include "trenchR.h"

double trenchDet(double *r,int n,double *v)
{
	int i;
	double logDet;
	
	logDet = 0.0;
	
	for (i = 1; i < n; i++)
		logDet += log(v[i]);

	return logDet;
}
	

