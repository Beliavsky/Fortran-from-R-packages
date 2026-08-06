// Matrix operations

#include "trenchR.h"

// Computation of dot product
double dot(int n,double* u,double* v)
{
   double t = 0.0;
   int i;
   
   for (i=0; i < n; i++)
      t += u[i] * v[i];
   
   return t;
}

// Vector - Matrix multiplication
void vecmat(int nc,double *v, double **m,double *u)
{
	int i,j;
	double nr;

	nr = nc; 
   
	for (i=0; i < nc; i++)
	{
	  double t = 0.0;
	  for (j = 0; j < nr; ++j)
		 t += v[j] * m[j][i];
	  u[i] = t;
	}

	return;
}
