// Convert list of vectors representing a doubly symmetric matrix to a matrix.

#include "trenchR.h"

void fromWedgeStorage(int n,double **B)
{
	int i,j;

	for (i = 1; i <= n; i++)
	for (j = 1; j <= n; j++)
		for (i = 1; i <= j; i++)
		{
			if (i > (int)((n+1-(j - i))/2))
				B[i-1][j-1] = B[n-j][n-i];
		}

	for (i = 0; i < n; i++)
		for (j = 0; j < i; j++)
		{	
			B[i][j] = B[j][i];
		}
			
	return;
}
