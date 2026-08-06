/**************************************************************************************
DLar
input parameters:
	c - vector of autocorrelations 1, r_1, ..., r_n
	nC - length of c, ie. n+1
	EPSL - machine epsilon in double precision

	output parameters:
	phikj - ar parameters, lags 1,...,n. Lag n+1 zero, ignore.
	phi - partial autocorrelations, lags 1,...,n. Lag n+1 zero, ignore.
	v - residual variances, lags 0, 1,...,n
	
	ifault:  0, OK;  1, c is not p.d.    
***************************************************************************************/
  
#include "nrutil.h"
#include <math.h>

void DLar(double *c, double *phikj, double *phi, double *v, int *nC, 
		double *EPSL, int *fault)  
{
	int i,j,k,n;
	double sum, EPS;

	VECTOR phiki;
    
	*fault = 0;
	n = *nC; 
	EPS = *EPSL;
	if (n < 1) *fault = 1;
	phiki = Vector(n);
	
	v[0] = c[0];
	
	if (c[0] <= EPS) *fault = 1;
	phi[0] = c[1] / c[0];
	phiki[0] = phi[0];
	phikj[0] = phi[0];
	v[1] = v[0] * (1.0 - phi[0]*phi[0]);
	
	if (v[1] <= EPS) *fault = 1;
	
	for (k = 2; k < n; k++)
	{
		sum = 0.0;
		for (i = 1; i < k; i++)
			sum += phiki[i - 1] * c[k - i];	         
		phi[k - 1] = (c[k] - sum) / v[k - 1];
		for (j = 1; j< k; j++)
			phikj[j - 1] = phiki[j - 1] - phi[k - 1] * phiki[k - j - 1];
		phikj[k - 1] = phi[k - 1];
		for (j = 1; j<=k; j++)
			phiki[j - 1] = phikj[j - 1];
		v[k] = v[k - 1] * (1.0 - phi[k - 1]*phi[k - 1]);
		if (v[k] <= EPS) *fault = 1;
		}
	
    
	free_vector(phiki);   
  
	return;    
}
