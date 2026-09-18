/*
    kzsv.c
    kz sample variance
    Copyright (C) 2005  Brian D. Close

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
 
#include <R.h>
#include "kz.h"

/* adaptive curvature function */
static double adaptive(double d, double m)
{
	/* m = max of the difference metric; a constant input gives m == 0 and
	   d/m = 0/0 = NaN, and (int)floor(q*NaN) downstream is undefined
	   behavior (it only survived via the min-window clamp by accident).
	   m <= 0 means no structure anywhere: keep the full window. */
	if (m <= 0.0) return 1.0;
	return( 1 - (d/m) );
}

SEXP kzsv(SEXP kza_data, SEXP kz_data, SEXP window, SEXP minimum_window_length, SEXP iterations, SEXP tolerance)
{
    SEXP ans=R_NilValue;
    
    PROTECT(ans = R_kzsv(kza_data, kz_data, window, minimum_window_length, tolerance));
    UNPROTECT(1);
    return(ans);
}

SEXP R_kzsv(SEXP kza_data, SEXP kz_data, SEXP window, SEXP minimum_window_length, SEXP tolerance)
{
    int i, t;
    int q;
    int qh, qt;
    int min_window_length;
    long n;
    double m;
    SEXP d, dprime;
    SEXP sigma=R_NilValue;
	double var, sum, avg;
	double temp;
	int size;
    double eps;
    
    eps = REAL(tolerance)[0];
	q = INTEGER_VALUE(window);
	min_window_length = INTEGER_VALUE(minimum_window_length);
    n = LENGTH(kza_data);
    
    if (LENGTH(kza_data)==0) error("KZA not specified.");
    if (LENGTH(kz_data)==0) error("KZ not specified.");
    if (LENGTH(kza_data)!=LENGTH(kz_data)) error("KZA KZ data size mismatch.");
    
	PROTECT(d = allocVector(REALSXP, n));
	PROTECT(dprime = allocVector(REALSXP, n));
    R_differenced(kz_data, d, dprime, q);

    m = R_maximum(d);

    PROTECT(sigma = allocVector(REALSXP, n));
    for (t=0; t<n; t++) {
	if (fabs(REAL(dprime)[t]) < eps) { /* dprime[t] = 0 */
	    qh = (int) floor(q*adaptive(REAL(d)[t], m));
	    qt = (int) floor(q*adaptive(REAL(d)[t], m));
	} else if (REAL(dprime)[t] < 0) {
	    qh = q;
   	    qt = (int) floor(q*adaptive(REAL(d)[t], m));
    	} else {
	    qh = (int) floor(q*adaptive(REAL(d)[t], m));
	    qt = q;
	}
	qt = ((qt) < min_window_length) ? min_window_length : qt;
	qh = ((qh) < min_window_length) ? min_window_length : qh;
	
        /* check bounds */
       	qh = (qh > n-t-1) ? n-t-1 : qh; /* head past end of series */
        qt = (qt > t) ? t : qt;
        
        /* variance */
        sum=0.0;
        for (i=t-qt; i<=t+qh; i++) {
            sum += REAL(kza_data)[i];
        }
        size = qh+qt+1;
        avg = sum/size;
        
        var=0.0;
        for (i=t-qt; i<=t+qh; i++) {
            temp = (REAL(kza_data)[i] - avg);
            var += temp * temp;
        }
        var /= size - 1;
        REAL(sigma)[t] = var;
    }
    UNPROTECT(3);
    return(sigma);
}

