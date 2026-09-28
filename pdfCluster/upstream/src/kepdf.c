#include<stdio.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>

void c_kepdfN(double *x, double *h, double *xeval, int *nx, int *ndim, int *neval, double *f, double *p)
{

 int i, j, k, _nx=*nx, _neval=*neval, _ndim=*ndim;
 double zero=0.0, s, p1=0.0, u=0.0, cons= pow( 0.3989423, _ndim )/_nx;

 for(i=0; i<_neval; i++)
 {
  s=zero;

   for(k=0; k<_nx; k++)
   {
    p1=zero;

    for(j=0; j<_ndim; j++)
    {
     u  = ( x[k+j*_nx]-xeval[i+j*_neval] )*h[k+j*_nx]; 
     p1+=  u*u;
    }
    s += exp(-p1)*p[k];
   }
    *f++=s*cons; 
 }
}




void c_kepdft7(double *x, double *h, double *xeval, int *nx, int
*ndim, int *neval, double *f, double *p)
{

 register int i, j, k;
 int _nx=*nx, _neval=*neval, _ndim=*ndim;
 double zero=0.0, one=1.0, seven=7.0, s, prod1, prod, u, cons=pow(0.3849915, _ndim)/_nx, seven_raised;  seven_raised = pow(seven, (double)_ndim*4.0);

 for(i=0; i<_neval; i++)
 {
  s=zero;

   for(k=0; k<_nx; k++)
   {
	prod = one;
	for(j=0; j<_ndim; j++)
{
	u  = ( *(x + k + j*_nx) - *(xeval + i + j*_neval) )**(h + k + j*_nx);
    prod*=(seven+u*u);
}
	prod1 = prod*prod;
	s += (seven_raised/(prod1*prod1))**(p+k);
   }
    *f++=s*cons;
 }
}


