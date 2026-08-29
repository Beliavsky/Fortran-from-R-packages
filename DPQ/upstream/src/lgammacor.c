/* Apart from these first lines: DIRECT COPY  from R sources  of
 *  <R>/src/nmath/lgammacor.c
 *      --------------------- @MM = ~/R/D/r-devel/R/src/nmath/lgammacor.c
 *
 *  Modified to have (nalgm, xbig) as arguments, plus .Call()able function:
 *
 *  Copyright (C) 2021 Martin Maechler
 */
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000-2021 The R Core Team
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
 *    #include <Rmath.h>
 *    double lgammacor(double x);
 *
 *  DESCRIPTION
 *
 *    Compute the log gamma correction factor ((for x >= 10)) so that
 *
 *    log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)
 *
 *    [ lgammacor(x) is called	Del(x)	in other contexts (dcdflib, toms708.c)]
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    written by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *  SEE ALSO
 *
 *    Loader(1999)'s stirlerr() {in ./stirlerr.c} is *very* similar in spirit,
 *    is faster and cleaner, but is only defined "fast" for half integers.
 */

#include "DPQpkg.h"

double dpq_lgammacor(double x, int nalgm, double xbig)
{
    const static double algmcs[15] = {  // R's  lgammacor() has nalgm = 5  hardwired
	+.1666389480451863247205729650822e+0,
	-.1384948176067563840732986059135e-4,
	+.9810825646924729426157171547487e-8,
	-.1809129475572494194263306266719e-10,
	+.6221098041892605227126015543416e-13,
	-.3399615005417721944303330599666e-15,
	+.2683181998482698748957538846666e-17,
	-.2868042435334643284144622399999e-19,
	+.3962837061046434803679306666666e-21,
	-.6831888753985766870111999999999e-23,
	+.1429227355942498147573333333333e-24,
	-.3547598158101070547199999999999e-26,
	+.1025680058010470912000000000000e-27,
	-.3401102254316748799999999999999e-29,
	+.1276642195630062933333333333333e-30
    };

/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 *   xbig = 2 ^ 26.5
 *   xmax = DBL_MAX / 48 =  2^1020 / 3 */

// R has these hardwired: ----------------
/* #define nalgm 5 */
/* #define xbig  94906265.62425156 */
#define xmax  3.745194030963158e306

/* In 'DPQ' pkg:  disable check to allow more experimentation :
    if (x < 10)
      /.******./
	ML_WARN_return_NAN
    else */ if (x >= xmax) {
	ML_WARNING(ME_UNDERFLOW, "lgammacor");
	/* allow to underflow below */
    }
    else if (x < xbig) { //  x in [10, xbig)  <==> t := 10/x  in [1.05e-7, 1]  <==> u:= 2t^2-1  in (-1, 1]
	double t = 10 / x;
	return chebyshev_eval(t * t * 2 - 1, algmcs, nalgm) / x;
    }
    return 1 / (x * 12);
}


// To be  .Call()ed  from R :
SEXP R_lgammacor(SEXP x_, SEXP nalgm_, SEXP xbig_)
{
    int nalgm = asInteger(nalgm_);
    if(nalgm <= 0)
	error("nalgm = %d <= 0", nalgm);
    else if(nalgm > 15)
	error("nalgm = %d > 15", nalgm);
    // else 1 <= nalgm <= 15 :
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    // other args are checked via  asReal(.), asInteger() below
    R_xlen_t i,	n = XLENGTH(x_);
    SEXP r_ = PROTECT(allocVector(REALSXP, n));
    double *x = REAL(x_), *r = REAL(r_),
	xbig = asReal(xbig_);
    for(i=0; i < n; i++) {
	r[i] = dpq_lgammacor(x[i], nalgm, xbig);
    }
    UNPROTECT(2);
    return r_;
}
