// Computing the inverse of a positive-definite Toeplitz matrix
// using the Trench algorithm

#include "trenchR.h"

int trenchInv(double *r,int n,double **b,double *v,double EPS)  
{
	double sigsq;
	int fault; 
	VECTOR fi;
	
	fi = Vector(n-1);
	fault = durlev (r,n,fi,v,&sigsq,EPS);
	
	if (fault != 0)
	{
		return fault;
		free_vector(fi);
	}

	v = Vector(n);
	upperWedge(fi,sigsq,n,b);
	
	free_vector(fi);

	return 0;
}








