/*************************************************************************************
trenchQFR

   Computing the Quadratic form and Determinant

   Takes arguments r,z, and EPS 
   Returns Quadraric Form S=y[0], Natural logarithm of determinant Log(Dn) = y[1],
   and an fault that describes the exit condition of trenchR.
   Possible values of fault and corresponding fault conditions are:
   0 The program is normally performed
   1 Error ("Singular Matrix")
   2 Error ("Input r[0] is not equal to 1.")
   3 Error ("The length(r) is not equal to the length(z))" 
***************************************************************************************/

#include "trenchR.h"

void trenchQFR(double *r,int *nn,
	 double *z,int *nnz,double *EPSL,double *y,int *fault)
{
	MATRIX b;
	VECTOR s,v;
	double logDet,_s,EPS;
	int i,n = *nn,nz = *nnz,_fault;
 
	if (n != nz)
	{
		for (i = 0; i < 2; i++)
			y[i] = 0.0;
		fault[0] = 3; 
		return;	
	}

	EPS = *EPSL;
	b = Matrix(n,n);
	v = Vector(n);

	_fault = trenchInv(r,n,b,v,EPS);

	if (_fault != 0)
	{
		for (i = 0; i < 2; i++)
			y[i] = 0.0;
		fault[0] = _fault; 
		free_matrix(b);
		free_vector(v);
		return;
	}
	else
		fault[0] = 0;

	fromWedgeStorage(n,b);
		
	logDet = trenchDet(r,n,v);

	s = Vector (n);

	vecmat(n,z,b,s);
	_s = dot(n,s,z);

	y[0] = _s;
	y[1] = logDet;
	
	free_matrix(b);
	free_vector(v);
	free_vector(s);
}

