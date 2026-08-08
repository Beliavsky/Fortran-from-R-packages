#ifndef BLASFUNS_H
#define BLASFUNS_H
#include <math.h>
static double cblas_ddot(int n,const double*x,int ix,const double*y,int iy){double s=0;for(int i=0;i<n;i++)s+=x[i*ix]*y[i*iy];return s;}
static void cblas_daxpy(int n,double a,const double*x,int ix,double*y,int iy){for(int i=0;i<n;i++)y[i*iy]+=a*x[i*ix];}
static void cblas_dscal(int n,double a,double*x,int ix){for(int i=0;i<n;i++)x[i*ix]*=a;}
static double cblas_dnrm2(int n,const double*x,int ix){return sqrt(cblas_ddot(n,x,ix,x,ix));}
static void cblas_dgemv(int order,int trans,int m,int n,double alpha,const double*A,int lda,const double*x,int ix,double beta,double*y,int iy){
  (void)order;(void)lda;
  if(trans==111){for(int i=0;i<m;i++){double s=0;for(int j=0;j<n;j++)s+=A[i*n+j]*x[j*ix];y[i*iy]=alpha*s+beta*y[i*iy];}}
  else {for(int j=0;j<n;j++){double s=0;for(int i=0;i<m;i++)s+=A[i*n+j]*x[i*ix];y[j*iy]=alpha*s+beta*y[j*iy];}}
}
#endif
