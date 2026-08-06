// Converts a doubly symmetric matrix to a list of vectors for more compact storage

#include "trenchR.h"
	
void upperWedge(double *phi,double sigsq,int n,double **b)
{
	int i,j,n1;
	VECTOR v = Vector(n);
	
	for (i = 1; i <= n-1; i++)
		v[i-1] = -phi[n-i-1] / sigsq;
	
	b[0][0] = 1.0 / sigsq;

	for (j = 2; j <= n; j++)
		b[0][j-1] = v[n-j];

	n1 = (int)((n-1)/2); 
	
	for (i = 2; i <= n1+1; i++)
		for (j = i; j <= n-i+1; j++)
			b[i-1][j-1] = b[i-2][j-2] + (v[n-j] * v[n-i] 
		   -v[i-2] * v[j-2]) * sigsq;
	
	free_vector(v);
	return;
}

