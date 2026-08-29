
/* --- this version, downloaded 2025-04-28  as  https://calgo.acm.org/1006.zip
**********************************************************************************

  Copyright (C) 2025 Martin Maechler --- for R | SEXP interface R_*() functions at end

  DELTAGAMMAINC Fast and Accurate Evaluation of a Generalized Incomplete Gamma
  Function.
  Copyright (C) 2016 Remy Abergel (remy.abergel AT gmail.com), Lionel Moisan (Lionel.Moisan AT parisdescartes.fr).

  This file is a part of the DELTAGAMMAINC software, dedicated to the
  computation of a generalized incomplete gammafunction. See the Companion paper
  for a complete description of the algorithm.

  ``Fast and accurate evaluation of a generalized incomplete gamma function''
  (Rémy Abergel, Lionel Moisan), preprint MAP5 nº2016-14, revision 1.

Then, "really published" as

  Rémy Abergel, Lionel Moisan. Algorithm 1006: Fast and accurate evaluation of a
  generalized incomplete gamma function. ACM Transactions on Mathematical Software, 2020, 46 (1), pp.1--24.
  DOI 10.1145/3365983. https://hal.science/hal-01329669v2

  This program is free software: you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************/

/* Definition of the G-function
 * ----------------------------
 *
 * We define the function G : (p,x) --> R as follows
 *
 * if x <= p:
 *
 * G(p,x) = exp(x-p*log(|x|)) * integral of s^{p-1} * exp(-sign(x)*s) ds from s = 0 to |x|
 *
 * otherwise:
 *
 * G(p,x) = exp(x-p*log(|x|)) * integral of s^{p-1} * exp(-s) ds from s = x to infinity
 *
 * where
 *
 * -  p is a positive real number
 * -  x is a real number, eventually equal to +infinity.
 *
 *
 * Definition of the generalized incomplete gamma function
 * -------------------------------------------------------
 *
 * We call generalized incomplete gamma function, and we note I_{x,y}^{mu,p},
 * the integral defined by
 *
 *   I_{x,y}^{mu,p} = integral of s^{p-1} * exp(-mu*s) ds
 *
 * where
 *
 * -  mu is a real number non equal to zero (in general we take mu = 1 or -1 but
 *  any nonzero real number is allowed)
 *
 * -  x and y are two nonnegative real numbers such as 0 <= x <= y <= +infinity,
 *    the setting y=+infinity is allowed only when mu > 0
 *
 * -  p is positive real number, p must be an integer when mu < 0.
 *
 *
 * * Brief description of several modules
 * --------------------------------------
 *
 *  -  Gammaln approximates the log of the complete gamma function
 *
 *      Gamma(p) = integral of s^{p-1} e^{-s} ds from s=0 to +infinity.
 *
 *     using Pugh's method (a refinement of Lanczos algorithm).
 *
 *  -  G_cfrac_lower evaluates the G-function in the domain x <= p using a
 *     continued fraction.
 *
 *  -  G_ibp evaluates the G-function in the domain (x < 0 and |x| < max(1,p-1))
 *     using a recursive integration by parts relation (WARNING: this function
 *     cannot be used when mu > 0).
 *
 *  -  G_cfrac_upper evaluates the G-function in the domain x > p using a
 *     continued fraction.
 *
 *
 * References
 * ----------
 *
 *-  R. Abergel and L. Moisan. 2016. Fast and accurate evaluation of a
 *   generalized incomplete gamma function, preprint MAP5 nº2016-14, revision 1
 *
 *-  F. W. J. Olver, D. W. Lozier, R. F. Boisvert, and C. W. Clark
 *   (Eds.). 2010. NIST Handbook of Mathematical Functions. Cambridge University
 *   Press. (see online version at [[http://dlmf.nist.gov/]])
 *
 *-  W. H. Press, S. A. Teukolsky, W. T. Vetterling, and
 *   B. P. Flannery. 1992. Numerical recipes in C: the art of scientific
 *   computing (2nd ed.).
 *
 *-  G. R. Pugh, 2004. An analysis of then Lanczos Gamma approximation (phd
 *   thesis)
 *
 */

#define DEBUG_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "DPQpkg.h" // -> "dpq.h"  and more; instead of R's  #include "nmath.h"

/* include "gamma_inc_T1006.h"  was kernel.h : -------------------
   MM FIXME: Some of these should be *optional* arguments  <<<<<<<<<<<<< FIXME ............
 */
/* several "macro variables" -- i.e. global tuning constants */

#define EPS 2.22e-16        /* Machine epsilon */
/* For the continued fractions G_cfrac_{lower,upper}() : */
#define ITMAX 1000000000    /* Maximum allowed number of iterations */
#define DPMIN 1.0e-300      /* Number near the smallest representable double-point number */

/* For Romberg integrations : */
#define NITERMAX_ROMBERG 15 /* Maximum allowed number of Romberg iterations */
#define TOL_ROMBERG 0.1     /* Tolerance factor used to stop the Romberg iterations */
#define TOL_DIFF 0.2        /* Tolerance factor used for the approximation of I_{x,y}^{mu,p} using differences */
/*--------------------------- end include gamma_inc_T... ------------------------*/

/* plim: compute plim(x), the frontier of the partition of the domain (p,x)
 * detailed in our paper.
 *
 *            |      x              if   0 < x
 *            |
 * plim(x) = <       0              if -9 <= x <= 0
 *            |
 *            | 5.*sqrt(|x|)-5.     otherwise
 * */
static double plim(double x)
{
  return (x >= 0.) ? x : ((x >= -9.) ? 0. : 5.*sqrt(-x)-5.);
}

/* *
 * Gammaln: compute log(Gamma(p)) for p > 0,
 *
 * where Gamma(p) = integral of s^{p-1} e^{-s} ds from s=0 to +infinity.
 *
 * The evaluation is done using the Pugh's method (approximation with 11 terms),
 * which is a refinement of the Lanczos method (approximation with 6 terms).
 *
 * input  : p is a positive real number.
 * output : log(Gamma(p))
 * */
double Gammaln(double p)
{
  double sum,z=p-1.;
  double d[11]={
    2.48574089138753565546e-5, 1.05142378581721974210e+0,
    -3.45687097222016235469e+0, 4.51227709466894823700e+0,
    -2.98285225323576655721e+0, 1.05639711577126713077e+0,
    -1.95428773191645869583e-1, 1.70970543404441224307e-2,
    -5.71926117404305781283e-4, 4.63399473359905636708e-6,
    -2.71994908488607703910e-9 };
  int k;
  for(sum=d[0],k=1; k<=10; k++)
      sum += d[k]/(z+(double)k);
  return (double)(log(1.860382734205265717*sum) -(z+0.5) + (z+0.5L)*log(z+11.400511));
}


/* *
 * G_cfrac_lower: compute G(p,x) in the domain x <= p using a continued fraction
 *
 * inputs :
 *
 *  -  p: positive real number
 *  -  x: a real number such as x <= p
 *
 * output : G(p,x)
 * */
double G_cfrac_lower(double p, double x)
{
  /* deal with special cases */
    if(x==0) return 0.;

  /* evaluate the continued fraction using Modified Lentz's method. However, as
   * detailed in our paper, we perform manually the first pass (n=1), of the
   * initial Modified Lentz's method. */
  double an=1., bn=p, // set an = a{1}, bn = b{1} (b{1} is non-zero)
      f= an/bn, c= an/DPMIN, d= 1./bn, del;
  int n=2;
  do {
      int k=n/2;
      an = ((n&1) ? k : -(p-1+k))*x; bn++;
      d=an*d + bn; if (d == 0) d=DPMIN;
      c=bn + an/c; if (c == 0) c=DPMIN;
      d=1.0/d;
      del=d*c;
      f *= del; n++;
  }
  while((fabs(del-1.0) >= EPS) && (n < ITMAX));

  return f;
}


/* *
 * G_ibp: compute G(p,x) in the domain (x<0) using a recursive integration by
 * parts
 *
 * inputs :
 *
 *  -  p : positive *integer*
 *  -  x : a negative real number such as |x| < max(1,p-1)
 *
 * output : G(p,x)
 * */
double G_ibp(double p, double x)
{
  double gln = Gammaln(p),
      t = fabs(x),
      tt = 1./(t*t);
  bool odd = (bool) ((int)p) % 2;

  /* main loop */
  double c = 1./t,
      d = p-1.,
      s = c*(t-d);
  int l = 0;
  bool stop;
  do {
    c *= d*(d-1)*tt;
    d -= 2.;
    double del = c*(t-d);
    s += del; l++;
    stop = (fabs(del) < fabs(s)*EPS);
  }
  while(!stop && (l < floor((p-2)/2)));

  if(odd && !stop) { s += d*c/t; }
  return (((odd)?-1.:1.) * exp(-t+gln-(p-1)*log(t)) + s)/t;
}


/* *
 * G_cfrac_upper: compute G(p,x) in the domain x > p using a continued fraction
 *
 * inputs :
 *
 *  -  p: positive real number
 *  -  x: a positive number x > p, eventually x=infinity
 *
 * output : G(p,x)
 * */
double G_cfrac_upper(double p, double x)
{
  /* deal with special cases */
    if(isinf(x)) return 0.;

  /* evaluate the continued fraction using Modified Lentz's method. However, as
   * detailed in our paper, we perform manually the first pass (n=1), of the
   * initial Modified Lentz's method. */
  double an=1., bn=x+1.0-p, // set an = a{1}, bn = b{1}
      c,d,del,f;
  bool t = (bn != 0);
  int n;
  if(t) { // b{1} is non-zero
    f=an/bn; c=an/DPMIN; d=1./bn;
    n=2;
  }
  else { // b{1}=0 but b{2} is non-zero, we compute Mcfrac = a{1}/f with f = a{2}/(b{2}+) a{3}/(b{3}+) ...
    an=-(1-p); bn=x+3.0-p;
    f=an/bn; c=an/DPMIN; d=1./bn;
    n=3;
  }
  int i = n-1;
  do {
    an = -i*(i-p); bn += 2.0;
    d=an*d + bn; if (d == 0) d=DPMIN;
    c=bn + an/c; if (c == 0) c=DPMIN;
    d=1.0/d; del=d*c; f*=del; i++; n++;
  }
  while((fabs(del-1.0) >= EPS) && (n < ITMAX));

  return (t) ? f : 1./f;
}


/* *
 * G_func: compute G(p,x) using the appropriate routine according to the value
 * of (p,x).
 *
 * inputs :
 *
 * -  p: a positive real number
 * -  x: a real number (eventually x=+inf)
 *
 * output : G(p,x)
 * */
double G_func(double p, double x)
{
    return (p >= plim(x)) ? G_cfrac_lower(p,x)
	// else p < plim(x) :
	: ((x < 0) ? G_ibp(p,x)
	/* x >= 0 */: G_cfrac_upper(p,x));
}


/* perform 1 iteration of the Romberg approximation of I_{x,y}^{mu,p} */
static void
romberg_iterations(double *R,
		   double sigma, int n, double x, double y, double mu, double p, double h, double pow2)
{
  int adr0_prev = ((n-1)*n)/2,
      adr0      = (n*(n+1))/2;

  double sum = 0.;
  for(int j=1; j <= pow2; j++) {
      double xx = x + ((y-x)*(2*j-1))/(2.*pow2);
      sum += exp(-mu*xx+(p-1)*log(xx)-sigma);
  }
  R[adr0] = 0.5*R[adr0_prev] + h*sum;

  /* main loop */
  double pow4 = 4.;
  for(int m=1; m<=n; m++) {
    R[adr0+m] = (pow4*R[adr0+(m-1)]-R[adr0_prev+(m-1)])/(pow4-1.);
    pow4 *= 4.;
  }

  return;
}

/* Estimate I_{x,y}^{mu,p} using a Romberg approximation: this algorithm
 * computes rho and sigma such as the Romberg approximation of I_{x,y}^{mu,p}
 * is given by I_{x,y}^{mu,p} = to rho * exp(sigma)                    */
void romberg_estimate(double *rho, double *sigma,
		      double x, double y, double mu, double p, double relerr_max, int niter_max)
{
  /* memory allocation */
  // was  R = (double*) malloc (((niter_max+1)*(niter_max+2))/2*sizeof(double));
  double *R = (double*) R_alloc(((niter_max+1)*(niter_max+2))/2,
				sizeof(double));
  /* initialization (n=1) */
  *sigma = -mu*y + (p-1)*log(y);
  R[0] = 0.5*(y-x)*(exp(-mu*x+(p-1)*log(x)-(*sigma))+1.);

  /* loop for n > 0 */
  int adr0 = 0,
      n = 1;
  double
      h = (y-x)/2., // n=1, h = (y-x)/2^n
      pow2 = 1. ,   // n=1; pow2 = 2^(n-1)
      relerr = ldexp(relerr_max, +1); // 2 * rel.max
  while ((n <= niter_max) && (relerr > relerr_max)) {
      // update R //
      romberg_iterations(R,*sigma,n,x,y,mu,p,h,pow2);
      h    /= 2;
      pow2 *= 2;
      // estimate relative error //
      adr0 = (n*(n+1))/2;  // ==> 2*niter_max < ~= sqrt(int_max)
      relerr = fabs((R[adr0+n]-R[adr0+n-1])/R[adr0+n]);
      n++;
  }

  /* save output (Romberg estimate) and free memory */
  *rho = R[adr0+(n-1)];
  /* free(R); */

  return;
}


/* deltagammainc: our algorithm for the approximation of I_{x,y}^{mu,p}.
 *
 * inputs :
 *
 * -  mu is a real number different from 0
 *
 * -  x and y are two nonnegative real numbers such as 0 <= x <= y <= +inf,
 *    and the setting y=+infinity is allowed only when mu > 0
 *
 * -  p is positive real number, p must be an integer when mu < 0
 *
 * outputs: this procedure computes (rho,sigma,method) described below:
 *
 * (rho,sigma) are such as the approximated value of I_{x,y}^{mu,p} is
 * I=(*rho)*exp(*sigma)
 *
 * (*method) is a flag describing how I was computed: (*method)==1 when I is
 * estimated using a difference, (*method)==2 when I is estimated using
 * Romberg's method, (*method) < 0 in other (trivial) cases (like I=0, or I =
 * Gamma(p)).
 * */
void deltagammainc(double *rho, double *sigma, int *method,
		   double x, double y, double mu, double p)
{
  /* deal with particular cases */
  if(x==y) {
      *rho = 0.;
      *sigma = -INFINITY;
      (*method) = -1;
      return;
  }
  if (x==0. && isinf(y)) {
      *rho = 1.;
      *sigma = Gammaln(p) - p*log(mu);
      (*method) = -2;
      return;
  }

  /* initialization */
  double
      mx = G_func(p, mu*x), nx = (isinf(x)) ? -INFINITY : -mu*x + p*log(x),
      my = G_func(p, mu*y), ny = (isinf(y)) ? -INFINITY : -mu*y + p*log(y);

  /* KERNEL: Compute (mA,nA) and (mB,nB) such as I_{x,y}^{mu,p} can be
     approximated by the difference A-B, where A >= B >= 0, A = mA*exp(nA) and B
     = mB*exp(nB). When the difference involves more than one digit loss due to
     cancellation errors, the integral I_{x,y}^{mu,p} is evaluated using the
     Romberg approximation method. */
  double mA,mB, nA,nB;
  if(mu < 0) {
      mA = my; nA = ny;
      mB = mx; nB = nx;
  } // mu >= 0
  else if(p < plim(mu*x)) {
      mA = mx; nA = nx;
      mB = my; nB = ny;
  }
  else if (p < plim(mu*y)) {
      double gln = Gammaln(p);
      mA = 1.; nA = gln-p*log(mu);
      nB = fmax(nx,ny);
      mB = mx*exp(nx-nB) + my*exp(ny-nB);
  }
  else {
      mA = my; nA = ny;
      mB = mx; nB = nx;
  }

  /* compute (rho,sigma) such as rho*exp(sigma) = A-B */
  *rho = mA - mB*exp(nB-nA);
  *sigma = nA;

  /* if the difference involved a significant loss of precision, compute Romberg estimate */
  if(!isinf(y) && ((*rho)/mA < TOL_DIFF)) {
      romberg_estimate(rho,sigma, x,y,mu,p,
		       /* relneeded = */ EPS/TOL_ROMBERG,
		       /* niter_max = */ NITERMAX_ROMBERG);
      (*method) = 2;
  }
  else (*method) = 1;
} // end{ deltagammainc }


// .Call(C_R_lgammaP11, x) :
SEXP R_lgammaP11(SEXP x_)
{
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    R_xlen_t n = XLENGTH(x_);
    SEXP r_ = PROTECT(allocVector(REALSXP, n));
    double *x = REAL(x_),
	*r = REAL(r_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = Gammaln(x[i]);
    UNPROTECT(2);
    return r_;
}


// .Call(C_R_deltagammainc, x,y, mu, p)  |-->  list(rho = *, sigma = *, method = *) :
SEXP R_deltagammainc(SEXP x_, SEXP y_, SEXP mu_, SEXP p_)
{
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    PROTECT(y_ = isReal(y_) ? y_ : coerceVector(y_, REALSXP));
    PROTECT(p_ = isReal(p_) ? p_ : coerceVector(p_, REALSXP));
    double mu = asReal(mu_);
    if(!R_FINITE(mu)) error("'mu' must be a finite numeric");

    R_xlen_t i, n,
	n_x = XLENGTH(x_),
	n_y = XLENGTH(y_),
	n_p = XLENGTH(p_);
    if(!n_x || !n_y || !n_p)
	n = 0;
    else { // otherwise, recycle to common length n :
	n = n_x;
	i = (n_y >= n_p) ? n_y : n_p;
	if(n < i) n = i;
    }
     // ==> n = max(length(x), length(y), length(p))

#ifdef DEBUG_
    /* if(trace) { */
    REprintf("R_deltagammainc(x,y, mu=%g, p): n := max(length(.), ..) = %lld\n", mu, (long long)n);
#endif

    SEXP rho_ = PROTECT(allocVector(REALSXP, n)),
	 sig_ = PROTECT(allocVector(REALSXP, n)),
	meth_ = PROTECT(allocVector(INTSXP , n));
    double *x = REAL(x_), *y = REAL(y_), *p = REAL(p_),
	*rho = REAL(rho_), *sig = REAL(sig_);
    int *meth = INTEGER(meth_);
    for(i=0; i < n; i++) {
	deltagammainc(rho+i, sig+i, meth+i,
		      x[i % n_x], y[i % n_y], mu,
		      p[i % n_p]);
    }
    const char *nms[] = {"rho", "sigma", "method", ""};
    SEXP ans = PROTECT(Rf_mkNamed(VECSXP, nms));
    SET_VECTOR_ELT(ans, 0, rho_);
    SET_VECTOR_ELT(ans, 1, sig_);
    SET_VECTOR_ELT(ans, 2, meth_);
    UNPROTECT(7);
    return ans;
}
