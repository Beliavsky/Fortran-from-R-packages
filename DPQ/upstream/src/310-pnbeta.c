/* 310-CumNonCentralBeta.f -- translated by f2c (version 20050501).
 * ----------------------- then, produced by
 * $Id: f2c-clean,v 1.10 2002/03/28 16:37:27 maechler Exp $
 */

#include "DPQpkg.h"

double betanc(double x, double a, double b, double lambda,
	      double errmax, int itrmax, int *ifault);

double ncbeta1(double a, double b, double lambda, double x,
	       Rboolean use_AS226,
	       int itrmax, /* only passed to betanc() for AS 226 currently,
			    *  was hardwired = 100 */
	       double errmax,
	       int *xj, // output: gives the number of iterations taken
	       int *ifault)
{
/*       ALGORITHM AS 310 APPL. STATIST. (1997), VOL. 46, NO. 1

       Computes the cumulative distribution function of a
       non-central beta random variable
*/

/*   Auxiliary Functions : ----------

 * <math.h> :

     double sqrt(double), log(double), exp(double), lgamma(double);

 * <Rmath.h> :

     double precision function betain(x, p, q, beta, ifault) --> pbeta()  ! (MM)
     double precision function gammad(x, p, ifault)          --> pgamma() ! (MM)

     double lgamma(a);
     double lbeta(a, b);
*/

    // Check for admissibility of parameters :

    if (lambda < 0. || a <= 0. || b <= 0.) { // MM: ncp = lambda = 0 is allowed !
	*ifault = 3; return x;
    }
    if (x < 0. || x > 1.) {
	*ifault = 2; return x;
    }
    if (x == 0. || x == 1.) {
	*ifault = 1; return x;
    }
    *ifault = 0;

    *xj = 0;

    // if (lambda < 54.) {/* AS 226 as it stands is sufficient in this situation */
    if (use_AS226) {

	return betanc(x, a, b, lambda, errmax, itrmax, ifault);

    } else { // !use_AS226 [was 'lambda >= 54' hardwired]

	double q, r__, s, s0, s1, t0, t1, fx, gx;
	double ebd, sum, lBeta, temp, psum;
	double errbd, ftemp;

	double c__ = lambda * 0.5;
	int m = (int) (c__ + 0.5);
	double t = 5. * sqrt((double) m);
	int
	    it = (int)ceil(t),
	    iterlo = m - it,
	    iterhi = m + it;

	t = -c__ + m * log(c__) - lgamma(m + 1.);
	q = exp(t);
	r__ = q;
	psum = q;
	// lBeta = lgamma(a + m) + lgamma(b) - lgamma(a + m + b);
	lBeta = lbeta(a + m, b); // from <Rmath.h>
	s1 = (a + m) * log(x) + b * log1p( - x) - log(a + m) - lBeta;
	gx = exp(s1);
	fx = gx;
	// r__1 = a + m;
/*	temp = betain_(x, &r__1, b, &beta, ifault); */
	// double pbeta(double x, double a, double b, int lower_tail, int log_p) :
	temp = pbeta(x, a + m, b, /*lower= */TRUE, /*log_p=*/FALSE);
	ftemp = temp;
	(*xj)++;
	sum = q * temp; // web version '310' had "SUM = Q - TEMP" which is incorrect.
	int iter1 = m;

/*      The first set of iterations starts from M and goes downwards */

	while(iter1 >= iterlo && q >= errmax) {
// L20:
	    q *= (iter1 / c__); //  '310', the online copy of AS 310 has "Q = Q - ITER1 / C" which is incorrect.
	    (*xj)++;
	    gx = (a + iter1) / (x * (a + b + iter1 - 1.)) * gx;
	    iter1--;
	    temp += gx;
	    psum += q;
	    sum += q * temp;
	}
// L30:
	t0 = lgamma(a + b) - lgamma(a + 1.) - lgamma(b);
	s0 = a * log(x) + b * log1p( - x);
	s = 0.; // << correction: initialization of 's' was missing in '310' *and* the publication (!)
	for (double j = 0.; j < iter1; ++j) {
	    s += exp(t0 + s0 + j * log(x));
	    t1 = log(a + b + j) - log(a + 1. + j) + t0;
	    t0 = t1;
	}

/*       Compute the first part of error bound */

	errbd = pgamma(c__, (double)iter1, 1.,/*lower= */FALSE, /*log_p=*/FALSE)
	    * (temp + s);
	q = r__;
	temp = ftemp;
	gx = fx;
	int iter2 = m;
	while(iter2 < iterhi && (ebd = errbd + (1. - psum) * temp) >= errmax) {
// L50:
	    iter2 ++;
	    (*xj)++;
	    q = q * c__ / iter2;
	    psum += q;
	    temp -= gx;
	    gx = x * (a + b + iter2 - 1.) / (a + iter2) * gx;
	    sum += q * temp;
	}
// L60:
	return sum;
    }
} /* ncbeta1 */

// call this via .C()  from R's pnbetaAS310() --> ../R/pnbeta.R
void ncbeta(double *a, double *b, double *lambda, double *x, int *n,
	    int *use_AS226,
	    double *errmax, int *itrmax,
	    int *ifault, double *res)
{
    int it_used;
    if(*ifault == 1) // all 4 main args are vectors of length n
	for(int i=0; i < *n; i++) // FIXME: *ifault depends on [i] -- must error() here
	    res[i] = ncbeta1(a[i], b[i], lambda[i], x[i], (Rboolean) use_AS226[i],
			     *itrmax, *errmax, &it_used, ifault);
    // TODO/FIXME : report 'it_used'  (for every 'i') if(verbose) ??

    else // x = x[1:n], the other three args are scalar :
	for(int i=0; i < *n; i++)
	    res[i] = ncbeta1(a[0], b[0], lambda[0], x[i], (Rboolean) use_AS226[0],
			     *itrmax, *errmax, &it_used, ifault);
    // TODO/FIXME : report 'it_used'  (for every 'i') if(verbose) ??
}


double betanc(double x, double a, double b, double lambda,
	      double errmax, int itrmax, int *ifault)
{
/*     ALGORITHM AS226 APPL. STATIST. (1987) VOL. 36, NO. 2
     Incorporates modification AS R84 from AS vol. 39, pp311-2, 1990

     Returns the cumulative probability of X for the non-central beta
     distribution with parameters A, B and non-centrality LAMBDA

     Auxiliary routines required: LGAMMA - log-gamma function (ACM
     291 or AS 245), and BETAIN - incomplete-beta function (AS 63)
*/
    /* Builtin functions */
    double sqrt(double), log(double), exp(double);

    /* Local variables */
    static double c__, q, a0, x0, ax, gx, lBeta, temp, sumq, errbd;

    /* Originally

    static real errmax = 1e-6f;
    static integer itrmax = 100;

    * Change ERRMAX and ITRMAX if desired ... */

    if (lambda < 0. || a <= 0. || b <= 0.) {
	*ifault = 2; return x;
    }
    if (x < 0. || x > 1.) {
	*ifault = 3; return x;
    }
    *ifault = 0;
    if (x == 0. || x == 1.) {
	return x;
    }

    c__ = lambda * 0.5;

    x0 = floor(fmax2(0., c__ - 5. * sqrt(c__)));
    a0 = a + x0;
    // lBeta = lgamma(a + x0) + lgamma(b) - lgamma(a0 + b);
    lBeta = lbeta(a + x0, b); // using <Rmath.h>  lbeta() which is better than the above
    temp = pbeta(x, a0, b, /*lower= */TRUE, /*log_p=*/FALSE);

    gx = exp(a0 * log(x) + b * log1p( - x) - lBeta - log(a0));
    if (a0 > a) {
	q = exp(-c__ + x0 * log(c__)) - lgamma(x0 + 1.);
    } else {
	q = exp(-c__);
    }
    int j = 0;
    ax = q * temp;
    sumq = 1. - q;
    double ret_val = ax;

/*     Recur over subsequent terms until convergence is achieved... */
// L10:
    do {
	j += 1;
	temp -= gx;
	gx = x * (a + b + j - 1.) * gx / (a + j);
	q = q * c__ / j;
	sumq -= q;
	ax = temp * q;
	ret_val += ax;

/*     Check for convergence and act accordingly... */
	errbd = (temp - gx) * sumq;
    } while (j < itrmax && errbd > errmax);

    if (errbd > errmax) {
	*ifault = 1;
    }
    return ret_val;
} /* betanc */

