// Durbin - Levinson algorithm (Golub, G. and Van Loan (1983). Matrix Computations, 
// John Hoptkins University Press, Baltimore, Algorithm 5.7-3.

#include "trenchR.h"

int durlev (double *c,int n,double *fi,double *v,double *vk,double EPS)
{
	int i, k, j;
	double s; 
	MATRIX phi;

	if (fabs(c[0] - 1.0) > EPS) // error("r[0] is not exactly equal to 1");
		return 2;
	
	phi = Matrix(n,n);

	phi[1][1] = c[1];
	v[0] = c[0];
	v[1] = 1.0 - phi[1][1] * phi[1][1];
	
	if (v[1] < EPS) // nrerror("Singular Matrix-1");
	{
		free_matrix(phi);
		return 1;
	}

	for (i = 2; i < n; i++)
	{
			for (k = i; k >= 1; k--)
				if (k == i)
				{
					s = 0.0;
					for (j = 1; j <= i - 1; j++)
						s += phi[i - 1][j] * c[i - j];
					phi[i][i] = (c[i] - s) * (1.0 / v[i - 1]);
					v[i] = v[i - 1] * (1.0 - phi[i][i] * phi[i][i]);
					if (v[i] < EPS) // error("Singular Matrix-2");
                    {
						free_matrix(phi);
						return 1;
					}
				}
				else
				{
					phi[i][k] = phi [i - 1][k] - phi[i][i] * phi[i - 1][i - k];
				}
	}
	

	for (i = 1; i < n; i++)
		fi[i-1] = phi[n-1][i];
	
	*vk = v[n-1];

	free_matrix(phi);

	return 0;
}
