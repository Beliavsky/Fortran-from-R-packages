/*  == expint: Exponential Integral and Incomplete Gamma Function ==
 *
 *  Functions to compute the incomplete gamma function
 *
 *     G(a,x) = int_x^infty t^{a-1} exp(-t) dt
 *
 *  for 'a' real and 'x' >= 0. [This differs from 'pgamma' of base R
 *  in that negative values of 'a' are admitted.]
 *
 *  Copyright (C) 2016-2026 Vincent Goulet
 *
 *  The code in part IMPLEMENTATION is derived from the GNU Scientific
 *  Library (GSL) v2.2.1 <https://www.gnu.org/software/gsl/>
 *
 *  Copyright (C) 2007 Brian Gough
 *  Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002 Gerard Jungman
 *
 *  The code in part R TO C INTERFACE is derived from R source code.
 *
 *  Copyright (C) 1998--2025 The R Core Team
 *  Copyright (C) 2003--2023 The R Foundation
 *  Copyright (C) 1995--1997 Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301, USA.
 *
 *  AUTHOR for the GSL: G. Jungman
 *  AUTHOR for expint: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *                     with much indirect help from the R Core Team
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "locale.h"
#include "expint.h"


/*
 *  IMPLEMENTATION OF THE WORKHORSE
 *
 *  Adapted from "special functions" material in the GSL.
 *
 */

/* Continued fraction which occurs in evaluation
 * of Q(a,x) or Gamma(a,x).
 *
 *              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
 *   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
 *             1 +   1 +     1 +   1 +      1 +   1 +
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
 *
 * Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
 * See gamma_inc_Q_CF() below.
 *
 */
double gamma_inc_F_CF(double a, double x)
{
    const int    nmax  =  5000;
    const double small =  R_pow_di(DBL_EPSILON, 3);

    double hn = 1.0;           /* convergent */
    double Cn = 1.0 / small;
    double Dn = 1.0;
    int n;

    /* n == 1 has a_1, b_1, b_0 independent of a,x,
       so that has been done by hand                */
    for (n = 2 ; n < nmax ; n++)
    {
	double an;
	double delta;

	if (E1_IS_ODD(n))
	    an = 0.5 * (n - 1)/x;
	else
	    an = (0.5 * n - a)/x;

	Dn = 1.0 + an * Dn;
	if (fabs(Dn) < small)
	    Dn = small;
	Cn = 1.0 + an/Cn;
	if (fabs(Cn) < small)
	    Cn = small;
	Dn = 1.0/Dn;
	delta = Cn * Dn;
	hn *= delta;
	if (fabs(delta-1.0) < DBL_EPSILON)
	    break;
    }

    if (n == nmax)
	warning(_("maximum number of iterations reached in gamma_inc_F_CF"));

    return hn;
}

/* Adapted from specfun/gamma_inc.c in GSL sources. Note that base R
 * function 'gammafn' and 'pgamma' are used for positive values of
 * 'a'. */
double gamma_inc(double a, double x)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(a))
	return a + x;
#endif

    if (x < 0.0)
	return(R_NaN);
    else if (x == 0.0)
	return gammafn(a);
    else if (a == 0.0)
	return expint_E1(x, 0);
    else if (a > 0.0)
	return gammafn(a) * pgamma(x, a, 1, 0, 0);
    else if (x > 0.25)
    {
	/* continued fraction seems to fail for x too small; otherwise
	   it is ok, independent of the value of |x/a|, because of the
	   non-oscillation in the expansion, i.e. the CF is
	   un-conditionally convergent for a < 0 and x > 0
	*/
	return exp((a - 1) * log(x) - x) * gamma_inc_F_CF(a, x);
    }
    else if (fabs(a) < 0.5)
    {
	/* expint: use the recursion for -0.5 < a < 0 (instead of a
	 * series expansion as in GSL), relying on the accuracy of
	 * pgamma for small values of 'a', but nevertheless treat
	 * this case separately to avoid rounding errors in the loop
	 * below */
	const double da = a + 1.0;
	const double gax = gammafn(da) * pgamma(x, da, 1, 0, 0);
	const double shift = exp(-x + a * log(x));

	return (gax - shift)/a;
    }
    else
    {
	/* a = fa + da; da >= 0 */
	const double fa = floor(a);
	const double da = a - fa;

	double gax  = (da > 0.0 ? gammafn(da) * pgamma(x, da, 1, 0, 0)
		                : expint_E1(x, 0));
	double alpha = da;

	/* Gamma(alpha-1,x) = 1/(alpha-1) (Gamma(a,x) - x^(alpha-1) e^-x) */
	do
	{
	    const double shift = exp(-x + (alpha - 1.0) * log(x));
	    gax = (gax - shift)/(alpha - 1.0);
	    alpha -= 1.0;
	} while (alpha > a);

	return gax;
  }
}


/*
 *  R TO C INTERFACE
 *
 *  Adapted from src/main/arithmetic.c in R sources and from a similar
 *  scheme in package actuar. Main difference: everything is in this
 *  one file.
 *
 */
#define mod_iterate1(n1, n2, i1, i2)            \
        for (i = i1 = i2 = 0; i < n;            \
             i1 = (++i1 == n1) ? 0 : i1,        \
             i2 = (++i2 == n2) ? 0 : i2,        \
             ++i)

/* Function called by .External() */
SEXP expint_do_gammainc(SEXP args)
{
    SEXP sx, sa, sy;
    R_xlen_t i, ix, ia, n, nx, na;
    double ai, *a, xi, *x, *y;
    Rboolean naflag = FALSE;

    args = CDR(args);	       /* drop function name from arguments */

    if (!isNumeric(CAR(args)) || !isNumeric(CADR(args)))
        error(_("invalid arguments"));

    na = XLENGTH(CAR(args));
    nx = XLENGTH(CADR(args));
    if ((na == 0) || (nx == 0))
        return(allocVector(REALSXP, 0));

    n = (nx < na) ? na : nx;

    PROTECT(sa = coerceVector(CAR(args), REALSXP));
    PROTECT(sx = coerceVector(CADR(args), REALSXP));
    PROTECT(sy = allocVector(REALSXP, n));
    a = REAL(sa);
    x = REAL(sx);
    y = REAL(sy);

    mod_iterate1(na, nx, ia, ix)
    {
        ai = a[ia];
        xi = x[ix];
	if (ISNA(ai) || ISNA(xi))
	    y[i] = NA_REAL;
        else if (ISNAN(ai) || ISNAN(xi))
	    y[i] = R_NaN;
        else
        {
	    y[i] = gamma_inc(ai, xi);
	    if (ISNAN(y[i])) naflag = TRUE;
        }
    }

    if (naflag)
        warning(R_MSG_NA);

    if (n == na)
        SHALLOW_DUPLICATE_ATTRIB(sy, sa);
    else if (n == nx)
        SHALLOW_DUPLICATE_ATTRIB(sy, sx);

    UNPROTECT(3);

    return sy;
}
