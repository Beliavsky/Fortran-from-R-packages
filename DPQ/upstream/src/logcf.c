#include "DPQpkg.h"

/*---------------------------------------------vvvvvvvvvvvvvvvvvvvvv-------
 * Cut'n'paste from R's sources  ~/R/D/r-devel/R/src/nmath/pgamma.c
 * [as of 2021-04-19, lines 67--118]
 * ADDING 'eps', 'maxit' trace option
 * using LDOUBLE where appropriate (such that eps < 1e-14 may work ??), 2025-05
 *-------------------------------------------------------------------------*/


/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
#define log2scalef 256
static // gmp::as.bigz(2)^256  (double prec = ldexp(1., 256); (==> -Wpedantic warning)
const double scalefactor = 115792089237316195423570985008687907853269984665640564039457584007913129639936.;


/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxiliary in log1pmx() (and hence bd0()) and lgamma1p()
 */
static
double logcf (double x, double i, double d, double eps /* ~ relative tolerance */,
	      int maxit, int trace)
{
    LDOUBLE c1 = (LDOUBLE)2 * d;
    LDOUBLE c2 = (LDOUBLE)i + d;
    LDOUBLE c4 = c2 + d;
    LDOUBLE a1 = c2;
    LDOUBLE b1 = i * (c2 - i * x);
    LDOUBLE b2 = (LDOUBLE)d * d * x;
    LDOUBLE a2 = c4 * c2 - b2;

#if 0
    assert (i > 0);
    assert (d >= 0);
#endif

    b2 = c4 * b1 - i * b2;

    int it = 0;
    while (FABS(a2 * b1 - a1 * b2) > FABS(eps * b1 * b2)) {
	LDOUBLE c3 = c2*c2*x;
	c2 += d;
	c4 += d;
	a1 = c4 * a2 - c3 * a1;
	b1 = c4 * b2 - c3 * b1;

	c3 = c1 * c1 * x;
	c1 += d;
	c4 += d;
	a2 = c4 * a1 - c3 * a2;
	b2 = c4 * b1 - c3 * b2;

        if(trace) REprintf("it=%2d: ==> |b2|=%g", it, fabs((double)b2));
	if(fabs((double)b2) > scalefactor) {
            if(trace) REprintf("  Lrg |b2|"); /* a1 /= scalefactor; .. etc */
	    a1 = LDEXP(a1, -log2scalef);
	    b1 = LDEXP(b1, -log2scalef);
	    a2 = LDEXP(a2, -log2scalef);
	    b2 = LDEXP(b2, -log2scalef);
	} else if(fabs((double)b2) < 1 / scalefactor) {
            if(trace) REprintf("  Sml |b2|"); /* a1 *= scalefactor;  .. etc */
	    a1 = LDEXP(a1, log2scalef);
	    b1 = LDEXP(b1, log2scalef);
	    a2 = LDEXP(a2, log2scalef);
	    b2 = LDEXP(b2, log2scalef);
	}
        if(trace) REprintf("|-> a2/b2=%.16g\n", (double)(a2 / b2));
        if(++it > maxit) {
            warning("non-convergence in logcf(), it = %d > maxit = %d", it, maxit);
	    break;
        }
    }
    if(trace && it <= maxit)
	REprintf("  logcf(*) used %d iterations.\n", it);
    return (double)(a2 / b2);
}
#undef log2scalef


// To be  .Call()ed  from R :
SEXP R_logcf(SEXP x_, SEXP i_, SEXP d_, SEXP eps_, SEXP maxit_, SEXP trace_)
{
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    // other args are checked via  asReal(.), asInteger() below
    R_xlen_t i,	n = XLENGTH(x_);
    SEXP r_ = PROTECT(allocVector(REALSXP, n));
    double *x = REAL(x_), *r = REAL(r_),
	ii = asReal(i_),
	d  = asReal(d_),
	eps= asReal(eps_);
    int maxit = asInteger(maxit_),
	trace = asInteger(trace_);
    if(ii <= 0.) error("i = %g <= 0", ii);
    if( d <  0.) error("d = %g <  0", d);

    for(i=0; i < n; i++) {
	r[i] = logcf(x[i], ii, d, eps, maxit, trace);
    }
    UNPROTECT(2);
    return r_;
}
