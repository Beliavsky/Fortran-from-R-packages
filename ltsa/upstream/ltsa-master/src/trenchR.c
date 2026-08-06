/**************************************************************************************
trenchR
   Takes arguments r and EPS
   Returns inverse matrix x and an fault that describes the exit condition of trenchR.
   Possible values of fault and corresponding fault conditions are:
	   0 The program is normally performed
	   1 Error ("Singular Matrix")
	   2 Error ("Input r[0] is not equal to 1") 
***************************************************************************************/
  
#include "trenchR.h"
void trenchR(double *r,int *nn,double *EPSL,
	double *x,int *fault)  
{
	MATRIX b;
	VECTOR v;
	int n,i,j,_fault;
	double EPS;
	
	n = *nn;
	EPS = *EPSL;

	b = Matrix(n,n);
        v = Vector(n);
	
	_fault = trenchInv(r,n,b,v,EPS);
    
	if (_fault != 0)
	{
		fault[0] = _fault; 
		free_matrix(b);
        free_vector(v);
		return;
	}
	else
		fault[0] = 0;

			
	fromWedgeStorage(n,b);
		
	for (i = 0;i < n;i++)
		for (j = 0; j < n; j++)
			x[i * n + j] = b[i][j];
		
    free_matrix(b);
    free_vector(v);
    return;    
}
