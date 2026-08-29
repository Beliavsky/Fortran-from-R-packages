/* Copyright (C) 2022-2023  Martin Maechler

   Apart from these first lines: DIRECT COPY  from R sources  of
   <R>/src/nmath/chebyshev.c
       --------------------- @MM = ~/R/D/r-devel/R/src/nmath/chebyshev.c  */
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2022 The R Core Team
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    int chebyshev_init(double *dos, int nos, double eta)
 *    double chebyshev_eval(double x, double *a, int n)
 *
 *  DESCRIPTION
 *
 *    "chebyshev_init" determines the number of terms for the
 *    double precision orthogonal series "dos" needed to insure
 *    the error is no larger than "eta".  Ordinarily eta will be
 *    chosen to be one-tenth machine precision.
 *
 *    "chebyshev_eval" evaluates the n-term Chebyshev series
 *    "a" at "x".
 *
 *  NOTES
 *
 *    These routines are translations into C of Fortran routines
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    Based on the Fortran routine dcsevl by W. Fullerton.
 *    Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).
 */

#include "DPQpkg.h"

/* NaNs propagated correctly */

/* returns typically  nos-1  unless has trailing zero coefficients : */
int chebyshev_init(const double dos[], int nos, double eta)
{
    double err = 0.;
    for (int i = nos-1; i >= 0; i--) {
	err += fabs(dos[i]);
	if (err > eta)
	    return i;
    }
    return -1;
}
/* NOTE: The above returns index  'n' such that  dos[0..n] should be used,
   ----  i.e., strictly, there are  n+1  terms !
*/


/* NOTE:  R's Mathlib chebyshev_eval() "gives"  T_0(x) == 1/2  instead of == 1,
   ----   or put differently,
          it computes      a_0/2 + a_1*T_1(x) + a_2*T_2(x) + .... + a_n*T_n(x)
*/

// Chebyshev Polynomial P(x; a[])
double chebyshev_eval(double x,         // typically in [-1, 1]
		      const double a[], // n coefficients a[0..(n-1)];  n <= length(a)
		      const int n)      // n-1 = (maximal) degree
{
#ifndef _not_in_DPQ_ // disable these to allow for more experimentation, incl. *extra*polation:
    if (n < 1    || n > 1000) ML_WARN_return_NAN;
    if (x < -1.1 || x > 1.1)  ML_WARN_return_NAN;
#endif
    double
	twox = ldexp(x, 1), //= x * 2
	b2 = 0, b1 = 0, b0 = 0;
    /* This is Clenshaw's algorithm, e.g. see Wikipedia, apart from just using  a[0]/2 as constant term  */
    for (int i = 1; i <= n; i++) {
	b2 = b1;
	b1 = b0;
	b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}


// To be  .Call()ed  from R :
SEXP R_chebyshev_eval(SEXP x_, SEXP a_, SEXP n_)
{
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    PROTECT(a_ = isReal(a_) ? a_ : coerceVector(a_, REALSXP));
    // n_ is checked via asInteger() below
    R_xlen_t i,	nx = XLENGTH(x_);
    SEXP r_ = PROTECT(allocVector(REALSXP, nx));
    double *x = REAL(x_), *a = REAL(a_), *r = REAL(r_);
    int n_a = asInteger(n_); /* <==> use a[0 .. n_a]  as coefficients */
    // n_a = -1 is allowed <==> empty polynomial === 0 {constant 0}
    n_a++; // chebyshev_eval(*, a, n_a)  will use  a[0..(n_a-1)]
    if(n_a > LENGTH(a_))
	error("n_a = %d > length(a) = %d", n_a, LENGTH(a_));
    for(i=0; i < nx; i++) {
	r[i] = chebyshev_eval(x[i], a, n_a);
    }
    UNPROTECT(3);
    return r_;
}


/** Determine the number of terms needed for a given coef[] vector
 ** @param eta  Above  "Ordinarily eta will be chosen to be one-tenth machine precision"
 **             R Mathlib's  nmath/<foo>.c  uses  eta = DBL_EPSILON/20  in both cases
 **/
SEXP R_chebyshev_nt(SEXP coef_, SEXP eta_)
{
    PROTECT(coef_ = isReal(coef_) ? coef_ : coerceVector(coef_, REALSXP));
    if(XLENGTH(coef_) > INT_MAX)
	error("length(%s) = %lld > max.int = %d", "coef", (long long)XLENGTH(coef_), INT_MAX);
    int n_a = chebyshev_init(REAL(coef_), LENGTH(coef_), asReal(eta_));
    UNPROTECT(1);
    return ScalarInteger(n_a);
}
