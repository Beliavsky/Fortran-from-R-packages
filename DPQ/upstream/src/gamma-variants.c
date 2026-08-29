/*
 * Originally a copy of <Rsrc>/src/nmath/gamma.c  (by R Core et al, see below)
 *                      ----------------------
 * The versions here provide *more* options, notably for experimentation,
 * notably on platforms such as the M1 mac where LDOUBLE is not really higher precision.
 *
 * These are
 *
 *  Copyright (C) 2022--2024 Martin Maechler,  maechler@stat.math.ethz.ch
 *
 * everything else, of course
 *
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000-2022 The R Core Team
 *  Copyright (C) 2002-2018 The R Foundation
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
 */

/*
 *  DESCRIPTION
 *
 *    This function computes the value of the gamma function.
 *
 *  NOTES
 *
 *    This function is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *    (e.g. http://www.netlib.org/slatec/fnlib/gamma.f)
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 *
 *    MM specialized the case of  n!  for n < 50 - for even better precision
 */

#include <float.h> /* DBL_MIN etc */

#include "DPQpkg.h" // -> "dpq.h"  and more; instead of R's  #include "nmath.h"

// For  lgammacor() -- our dpq_lgammacor()  --> ./lgammacor.c :
#define lgcor_nalgm 5
#define lgcor_xbig  94906265.62425156
#define lgammacor(_Y_)  dpq_lgammacor(_Y_, lgcor_nalgm, lgcor_xbig)


double gammafn_ver(double x, int version, int trace_lev, stirlerr_version_t stirl_ver)
/* version: integer >= 1  for gamma() version:

    -  R version        | svn r | Date       | comment
    --------------------|-------|------------|----------
    1: 1.0.0 -1.7.0     |  7982 | 2000-02-05 |
    2: 1.7.1(?)- 2.5.x  | 24344 | 2003-05-16 | gamma(n)
    3: 2.6.0-           | 42549 | 2007-08-18 | ML_UNDERFLOW
    4: 3.6.0 - 4.3.x    | 75839 | 2018-12-12 | no Inf/0 warning
    5: 4.4.0 ?          | ?     | 2024-..    | accuracy

================ UNFINISHED ===========================
 now contains mostly  the R-devel version '5'
================ UNFINISHED ===========================
*/
{
    const static double gamcs[42] = {
	+.8571195590989331421920062399942e-2,
	+.4415381324841006757191315771652e-2,
	+.5685043681599363378632664588789e-1,
	-.4219835396418560501012500186624e-2,
	+.1326808181212460220584006796352e-2,
	-.1893024529798880432523947023886e-3,
	+.3606925327441245256578082217225e-4,
	-.6056761904460864218485548290365e-5,
	+.1055829546302283344731823509093e-5,
	-.1811967365542384048291855891166e-6,
	+.3117724964715322277790254593169e-7,
	-.5354219639019687140874081024347e-8,
	+.9193275519859588946887786825940e-9,
	-.1577941280288339761767423273953e-9,
	+.2707980622934954543266540433089e-10,
	-.4646818653825730144081661058933e-11,
	+.7973350192007419656460767175359e-12,
	-.1368078209830916025799499172309e-12,
	+.2347319486563800657233471771688e-13,
	-.4027432614949066932766570534699e-14,
	+.6910051747372100912138336975257e-15,
	-.1185584500221992907052387126192e-15,
	+.2034148542496373955201026051932e-16,
	-.3490054341717405849274012949108e-17,
	+.5987993856485305567135051066026e-18,
	-.1027378057872228074490069778431e-18,
	+.1762702816060529824942759660748e-19,
	-.3024320653735306260958772112042e-20,
	+.5188914660218397839717833550506e-21,
	-.8902770842456576692449251601066e-22,
	+.1527474068493342602274596891306e-22,
	-.2620731256187362900257328332799e-23,
	+.4496464047830538670331046570666e-24,
	-.7714712731336877911703901525333e-25,
	+.1323635453126044036486572714666e-25,
	-.2270999412942928816702313813333e-26,
	+.3896418998003991449320816639999e-27,
	-.6685198115125953327792127999999e-28,
	+.1146998663140024384347613866666e-28,
	-.1967938586345134677295103999999e-29,
	+.3376448816585338090334890666666e-30,
	-.5793070335782135784625493333333e-31
    };

#ifdef NOMORE_FOR_THREADS
    static int ngam = 0;
    static double xmin = 0, xmax = 0., xsml = 0., dxrel = 0.;

    /* Initialize machine dependent constants, the first time gamma() is called.
       FIXME for threads ! */
    if (ngam == 0) {
	ngam = chebyshev_init(gamcs, 42, DBL_EPSILON/20);/*was .1*d1mach(3)*/
	gammalims(&xmin, &xmax);/*-> ./gammalims.c */
	xsml = exp(fmax2(log(DBL_MIN), -log(DBL_MAX)) + 0.01);
	/*   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE */
	dxrel = sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)) */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 * (xmin, xmax) are non-trivial, see ./gammalims.c
 * xsml = exp(.01)*DBL_MIN
 * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
*/
# define ngam 22
// now allowing denormalized result (2^-1074) for xmin *and* round to be slightly *outside*
/* # define xmin -182.     // -177.56341258681965  was -170.5674972726612 */
# define xminNorm -170.5674972726612
/* # define xmax  171.6244 //  171.62437695630271  was  171.61447887182298 */
# define xsml 2.2474362225598545e-308
/* # define dxrel 1.490116119384765696e-8 -- set below, version-dependently */
    // *variable* as they change with version :
    static double xmin = 0., xmax = 0., dxrel = 0.;
#endif

    if(ISNAN(x)) return x;

    if(trace_lev >= 1) REprintf("gamma_ver(%g, version=%d): ", x, version);

    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if(x == 0 || (x < 0 &&
		  ((version <= 2 && x == (long long)x) ||
		   /* later versions */ x == round(x)))) {
	ML_WARNING(ME_DOMAIN, "gammafn_ver");
	return ML_NAN;
    }

    if(version <= 4) { // Version from R 1.0.0 till 4.3.x (at least)
	xmin = -170.5674972726612;
	xmax =  171.61447887182298;
    } else { // version == 5 (latest for now)
	xmin = -182.;
	xmax =  171.6244;
    }

    if (x > xmax) {			/* Overflow */
	if(trace_lev >= 1) REprintf("x > xmax --> +Inf\n");
	if(version <= 3)
	    ML_WARNING(ME_RANGE, "gammafn_ver");
	// else No warning: +Inf is the best answer
	return ML_POSINF;
    }

    if (x < xmin) {			/* Underflow */
	if(trace_lev >= 1) REprintf("x < xmin --> 0.\n");
	if(version <= 3) { // ver <= 2: DBL_MIN*DBL_MIN; ver = 3: 0 + warning
	    ML_WARNING(ME_UNDERFLOW, "gammafn_ver");
	    return (version == 3) ? 0. : (DBL_MIN * DBL_MIN);
	    // had  #define ML_UNDERFLOW  (DBL_MIN * DBL_MIN)   in <nmath.h>
	}
	else // No warning: 0 is the best answer
	    return 0.;
    }

    // *only* influences  if  "full precision was not achieved"  warning happens
    if(version <= 2)
	dxrel = 67108864.; // dxrel = sqrt(1/DBL_EPSILON) = 2 ^ 26
    else /* r44772 | ripley | 2008-03-17 : ==> R 2.7.0
	  * dxrel was the reciprocal of the correct value */
	dxrel = 1.490116119384765696e-8; // = sqrt(DBL_EPSILON) = 2^{-26}

    if(trace_lev >= 1) REprintf("x 'normal'; dxrel = %g\n", dxrel);

    /* Compute gamma(x) for  xmin <= x <= xmax */

  if(version <= 4) {
    int i, n;
    double sinpiy, value;
    double y = fabs(x);

    if (y <= 10) {

	/* Compute gamma(x) for -10 <= x <= 10
	 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
	 * first of all. */

	n = (int) x;
	if(x < 0) --n;
	y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */

	if(trace_lev >= 1) REprintf("version <= 4, |x| <= 10: n=%d, y=%g\n", n, y);
	--n;
	value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
	if(trace_lev >= 1) REprintf("version <= 4, |x| <= 10: n=%d, y=%g --> val. := chebyshev(*)= %g\n",
				    n, y, value);
	if (n == 0)
	    return value;/* x = 1.dddd = 1+y */

	if (n < 0) {
	    /* compute gamma(x) for -10 <= x < 1 */

	    /* exact 0 or "-n" checked already above */

	    /* The answer is less than half precision */
	    /* because x too near a negative integer. */
	    if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
		ML_WARNING(ME_PRECISION, "gammafn_ver");
	    }

	    /* The argument is so close to 0 that the result would overflow. */
	    if (y < xsml) {
		ML_WARNING(ME_RANGE, "gammafn_ver");
		if(x > 0) return ML_POSINF;
		else return ML_NEGINF;
	    }

	    n = -n;

	    for (i = 0; i < n; i++) {
		value /= (x + i);
	    }
	    return value;
	}
	else {
	    /* gamma(x) for 2 <= x <= 10 */

	    for (i = 1; i <= n; i++) {
		value *= (y + i);
	    }
	    return value;
	}
    }
    else {
	/* gamma(x) for	 y = |x| > 10. */

	if(trace_lev >= 1) REprintf("version <= 4, |x| > 10 ..\n");

	if (x > xmax) {			/* Overflow */
	    // No warning: +Inf is the best answer
	    return ML_POSINF;
	}

	if (x < xmin) {			/* Underflow */
	    // No warning: 0 is the best answer
	    return 0.;
	}

	if(version >= 2 &&
	   y <= 50 && y == (int)y) { /* compute (n - 1)! */
	    value = 1.;
	    for (i = 2; i < y; i++) value *= i;
	}
	else { /* normal case */
	    value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI +
			((2*y == (int)2*y)? dpq_stirlerr(y, stirl_ver) : lgammacor(y)));
	}
	if (x > 0)
	    return value;

	if (fabs((x - (int)(x - 0.5))/x) < dxrel) {

	    /* The answer is less than half precision because */
	    /* the argument is too near a negative integer. */

	    ML_WARNING(ME_PRECISION, "gammafn_ver");
	}

	sinpiy = sinpi(y);
	if (sinpiy == 0) {		/* Negative integer arg - overflow */
	    ML_WARNING(ME_RANGE, "gammafn_ver");
	    return ML_POSINF;
	}

	return -M_PI / (y * sinpiy * value);
    }

  } else { // version == 5 (latest for now; was in my R-devel for a while; no longer; e.g. LDOUBLE "not there" on M1

	int i, n = (int) x;

	if(trace_lev >= 1) REprintf("version >= 5, n=%d\n", n);
	if(n == 1) { //  x = 1.dddd = 1+y    //  ldexp(y, 1) == y*2

	    if(trace_lev >= 1) REprintf("n=1 -> return chebyshev(.)\n");
	    return chebyshev_eval(ldexp((double)((LDOUBLE)x - n), 1) - 1, gamcs, ngam) + .9375;
	}
	if(x < 0) --n; // <==> n := floor(x)
	LDOUBLE x_ = (LDOUBLE) x;
	LDOUBLE y = x_ - n;                //  ==> y in [0, 1)
	--n; // = floor(x) - 1
	LDOUBLE value = chebyshev_eval(ldexp((double)y, 1) - 1, gamcs, ngam) + .9375;
	/*
	 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
	 * first of all. */
	if(trace_lev >= 1) REprintf("y=%" PR_g_ ", n=%d, value=%" PR_g_,
				    y, n, value);
	if (n < 0) {
	    /* compute gamma(x) for xmin <= x < 1 */

	    /* exact 0 or "-n" checked already above */

	    /* The answer is less than half precision */
	    /* because x too near a negative integer. */
	    if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
		ML_WARNING(ME_PRECISION, "gammafn_ver");
	    }

	    /* The argument is so close to 0 that the result would overflow. */
	    if (y < xsml) {
		ML_WARNING(ME_RANGE, "gammafn_ver");
		if(x > 0) return ML_POSINF;
		else return ML_NEGINF;
	    }

	    if(!(sizeof(LDOUBLE) > sizeof(double))) { // <==> !capabilities("long.double")
		if (x < xminNorm) { /* x \in [xmin, xminNorm) ~= [-180, -170]
				       i.e., the region where gamma(x) is subnormal > 0 */
		    // using the code from lgamma() :
		    double y = fabs(x), sinpiy = fabs(sinpi(y));
		    if(trace_lev >= 1) REprintf("n < 0, *not* long.double, x < xminNorm=%g: -> |sinpi(y)|=%g",
						xminNorm, sinpiy);
		    return /* sign : gamma(x) is negative if (x % 2) > 1 <==> n is even  */
			(n % 2 == 0 ? -1. : 1.) *
			exp(M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x
			    - log(sinpiy) - lgammacor(y));
		}
	    }
	    // else x >= xminNorm  or have accurate  LDOUBLE :
	    n = -n;

	    if(trace_lev >= 1) REprintf("n < 0, long.double --> value /= prod_{i=0}^{%d-1} (x_ - i)\n%s", n,
					(trace_lev >= 2)? "value = (" : "");
	    for (i = 0; i < n; i++) {
		value /= (x_ + i);
 		if(trace_lev >= 2) REprintf("%" PR_g_ "%s", value, (i < n-1) ? ", " : "");
	    };  if(trace_lev >= 2) REprintf(")\n");
	}
	else { // n >= 1  <==> floor(x) >= 2 :  gamma(x) for 2 <= x <= xmax */
	    if(trace_lev >= 1) REprintf("n >= 1 --> value *= prod_{i=1}^{%d} (x_ - i)\n%s", n,
					(trace_lev >= 2)? "value = (" : "");
	    for (i = 1; i <= n; i++) {
		value *= (x_ - i); // was  *= (y + i)
		if(trace_lev >= 2) REprintf("%" PR_g_ "%s", value, (i < n) ? ", " : "");
	    };  if(trace_lev >= 2) REprintf(")\n");

	}
	return (double) value;
    } // version >= 5
}


SEXP R_gamma_ver(SEXP x_, SEXP version_, SEXP trace_, SEXP stirlerr_v)
{
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    R_xlen_t i, n = XLENGTH(x_);
    int i_ver   = asInteger(version_),
	i_trace = asInteger(trace_);
    SEXP ans = PROTECT(allocVector(REALSXP, n));
    double *x = REAL(x_), *r = REAL(ans);
    stirlerr_version_t stirl_ver = (stirlerr_version_t) asInteger(stirlerr_v);

    for(i=0; i < n; i++) {
	r[i] = gammafn_ver(x[i], i_ver, i_trace, stirl_ver);
    }
    UNPROTECT(2);
    return ans;
}
