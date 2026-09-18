#include <R.h>
#include "kz.h"

double R_maximum(SEXP v) 
{
	double m;
	int i;
	
	for(i=0, m = REAL(v)[0]; i<LENGTH(v); i++) {if (R_FINITE(REAL(v)[i])) m = MAX(REAL(v)[i], m);}
	return m;
}

void R_differenced(SEXP y, SEXP d, SEXP dprime, int q)
{
    int i;
    long n;
    
    n = LENGTH(y);    

	/* calculate d = |Z(i+q) - Z(i-q)| */
	for (i=0; i<q; i++) {REAL(d)[i] = fabs(REAL(y)[i+q] - REAL(y)[0]);}
	for (i=q; i<n-q; i++) {REAL(d)[i] = fabs(REAL(y)[i+q] - REAL(y)[i-q]);}
	for (i=n-q; i<n; i++) {REAL(d)[i] = fabs(REAL(y)[n-1] - REAL(y)[i-q]);}

	/* d'(t) = d(i+1)-d(i) */
	for(i=0; i<n-1; i++) REAL(dprime)[i] = REAL(d)[i+1]-REAL(d)[i];
	REAL(dprime)[n-1] = REAL(dprime)[n-2];
}

void differenced(double y[], double d[], double dprime[], long n, int q)
{
    int i;
    
	/* calculate d = |Z(t+q) - Z(t-q)| */
	for (i=0; i<q; i++) {d[i] = fabs(y[i+q] - y[0]);}
	for (i=q; i<n-q; i++) {d[i] = fabs(y[i+q] - y[i-q]);}
	for (i=n-q; i<n; i++) {d[i] = fabs(y[n-1] - y[i-q]);}

	/* d'(t) = d(t+1)-d(t) */
	for(i=0; i<n-1; i++) dprime[i] = d[i+1]-d[i];
	dprime[n-1] = 0;
}
