/* -*- mode: C ; delete-old-versions: never -*-

 Original fortran code from preprint
 (then via citeseer, now ~/save/papers/Numerics/wienergerm-approximation.pdf )
                         ------------------------------------------------------
    Spiridon Penev and Tenko Raykov (1997).
    A Wiener Germ Approximation of the Noncentral Chi Square
    Distribution and of Its Quantiles

    Report No. S97-6, (July 1997);
    Department of Statistics, School of Mathematics,
    The University of New South Wales, Australia

 Then published (((unchanged, just got new page numbers,
                   as seen from scanned PDF available from Researchgate
		   == ~/save/papers/Numerics/wienergerm-approx-2000.pdf )))
 as

 @article{PenR2000,
   author =	 {Penev, Spiridon and Raykov, Tenko},
   year =	 2000,
   title =	 {A {W}iener Germ approximation of the noncentral chi square
		   distribution and of its quantiles},
   journal =	 {Computational Statistics},
   volume =	 15,
   number =	 2,
   pages =	 {219--228},
   issn =	 {0943-4062},
   doi =		 {10.1007/s001800000029},
   url =		 {http://dx.doi.org/10.1007/s001800000029},
   publisher =	 {Physica-Verlag},
   keywords =	 {Key words: Noncentral chi square, Wiener Germ Approximation},
   language =	 {English}
 }

Penev, Spiridon and Raykov, Tenko (2000).
Computational Statistics \bold{15}, 219--228. \doi{10.1007/s001800000029}

   translated by f2c (version 20031025) and by f2c-clean,v 1.10;
   and simplified by Martin Maechler, Jan.2004
*/

#include "DPQpkg.h" // before stdio (for MINGW_...)

#include <float.h> /* DBL_MIN etc */

#include <R_ext/Print.h>

/* R: use pnorm(z) instead of  derfc(- z/sqrt(2)) / 2  ! */

static double h(double y)
{
/* h(y) := ( (1-y)* log(1-y) + y - y^2/2) / y^2
 *  --
 * numerically stable; Martin Maechler (2004) */
    static const double // -pedantic warns / clang *error*s about sqrt(.) , pow(.) not being const.
	E_2 = 2.7205673729636575453329564e-8, // == sqrt((10./3.) * DBL_EPSILON),// ~= 2.72  e-8
	E_3 = 1.0354681459788919721672346e-5; // == pow(5. * DBL_EPSILON, 1./3.);// ~= 1.035 e-5

    if(y == 1.) return(0.5);

    double ay = fabs(y);

#ifdef DEBUG_h
    Rprintf("h(%g)\n", y);
#endif
    if(ay < E_2) return (y * (1. +  y/2.)/6.);
    if(ay < E_3) return (y * (1./6. + y*(1./12. + y/20.)));
    // else: |y| >= E3 ~= 1.035e-5
    return (((1.-y)*log1p(- y) + y*(1.- y/2.)) / (y*y));
}

double nonc_chi(double x, /* int *n, */ /* vectorized in x : x[1..n] : */
		double ncp, double df, logical lower_tail, logical log_p,
		int variant)
{
    double mu2, ff, m2s, d_h, s, s_1, eta, h_1s, z, f2, d__2;

/* check of initial values */

    if (x  <= 0.) { return R_DT_0; }

/* start calculation */

    mu2 = ncp / df;
    ff = sqrt(1. + 4. * x * mu2 / df); /* === 1. + 2 * mu2 * s */
    m2s = 0.5 * (ff - 1.); /* === mu2 * s */
    /* FIXME:	 ^^^^^^ use Taylor approx. for ff ~= 1
     *		(<==> small x, ncp or large df) */
    s = m2s / mu2;
    if (s == 1.) {// MM: variant 1 works!
	// ==> do something far better than the proposed
	// p = (log_p ? -M_LN2 : 0.5); return;
	z = s_1 = h_1s = 0.;
	if(variant != 1) {
	    REprintf("nonc_chi(*, variant=%d) -> s = 1 ==> setting variant := 1\n",
		     variant);
	    variant = 1;
	}
    } else {
	// FIXME: when s is very close to 1, the following loses lots of precision!
	s_1 = s - 1.;
	h_1s = -h(1. - 1. / s);/* == h(1 - s) ;
				* MM: FIXME: depends on s which form is better */
	z = df * s_1*s_1 * (0.5 / s + mu2 - h_1s / s)
	    - log((1. - 2. * h_1s / ff) / s);

    }

#ifdef DEBUG_h
    Rprintf("c(ff= %g, s= %g, h= %g, z=%g)\n", ff,s, h_1s, z);
#endif

    if (variant == 1) {/* variant = 'f' : improved first order approximation */
	d__2 = 1. + 3 * mu2;
	z += 2./9. * d__2 * d__2 / (df * pow(1. + 2 * mu2, 3.));
    }
    else { /* variant = 's' : the second order approximation */

	f2 = ff * ff;
	d__2 = 1. + 3 * m2s;
	d_h = -1.5 * (1. + 4 * m2s)/ f2 + 5./3.* d__2*d__2 /(ff*f2);

	eta = ff - 2. * h_1s;
	eta = (eta - s * ff) / eta;

	ff = eta / (s_1 * s_1) / ff;
	z += 2*(d_h + 2 * d__2 / s_1 / f2 +
		ff* (3. - (0.5 + h(eta)) * eta)) / df;
    }

#ifdef DEBUG_h
    Rprintf(" zz= %g)\n", z);
#endif
    z = sqrt(fabs(z));/* z < 0 is possible : e.g. x=ncp=2, df=1 */
    if (s < 1.)
	z = -z;

    return pnorm(z, 0., 1., lower_tail, log_p);
} /* nonc_chi */

// Called via .C():
void pchisqV(double *x, int *n, /* vectorized in x : x[1..n] : */
	     double *ncp, double *df,
	     logical *lower_tail, logical *log_p, int *variant)
{
    if (*ncp < 0.) { error("pchisqV(): ncp < 0 is invalid\n"); return; /*-Wall*/ }
    if (*df  < 0.) { error("pchisqV(): df < 0 is invalid\n"); return; /*-Wall*/ }
    /* df = 0, ncp > 0 is allowed, so,
       in principle also the \delta_0 point mass for df = ncp = 0 */

    for(int i = 0; i < *n; i++) {
	x[i] = nonc_chi(x[i], *ncp, *df, *lower_tail, *log_p, *variant);
    }
    return;
} /* pchisqV */
