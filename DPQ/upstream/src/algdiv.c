#include "DPQpkg.h"

/*---------------------------------------------vvvvvvvvvvvvvvvvvvvvv-------
 * Cut'n'paste from R's sources  ~/R/D/r-devel/R/src/nmath/toms708.c
 * [as of 2019-08-15, lines 2182--2244]
 * replacing s/alnrel/log1p/ in one place
 *-------------------------------------------------------------------------*/
double algdiv(double a, double b)
{
/* ----------------------------------------------------------------------- */

/*     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8 */

/*                         -------- */

/*     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY */
/*     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X). */

/* ----------------------------------------------------------------------- */

    /* Initialized data */

    static double c0 = .0833333333333333;
    static double c1 = -.00277777777760991;
    static double c2 = 7.9365066682539e-4;
    static double c3 = -5.9520293135187e-4;
    static double c4 = 8.37308034031215e-4;
    static double c5 = -.00165322962780713;

    double c, d, h, t, u, v, w, x, s3, s5, x2, s7, s9, s11;

/* ------------------------ */
    if (a > b) {
	h = b / a;
	c = 1. / (h + 1.);
	x = h / (h + 1.);
	d = a + (b - 0.5);
    }
    else {
	h = a / b;
	c = h / (h + 1.);
	x = 1. / (h + 1.);
	d = b + (a - 0.5);
    }

/* Set s<n> = (1 - x^n)/(1 - x) : */

    x2 = x * x;
    s3 = x + x2 + 1.;
    s5 = x + x2 * s3 + 1.;
    s7 = x + x2 * s5 + 1.;
    s9 = x + x2 * s7 + 1.;
    s11 = x + x2 * s9 + 1.;

/* w := Del(b) - Del(a + b) */

    t = 1./ (b * b);
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 *
	    s3) * t + c0;
    w *= c / b;

/*                    COMBINE THE RESULTS */

    u = d * log1p(a / b); // R's toms708.c uses its own  alnrel(.) instead
    v = a * (log(b) - 1.);
    if (u > v)
	return w - v - u;
    else
	return w - u - v;
} /* algdiv */
// ---------  ------------------------- end cut'n'paste from R sources

// To be  .Call()ed  from R :
SEXP R_algdiv(SEXP a_, SEXP b_)
{
    if(!isReal(a_) || !isReal(b_))
	error("'a' and 'b' must be \"double\" numeric vectors");
    R_xlen_t i,
	n_a = XLENGTH(a_),
	n_b = XLENGTH(b_),
	n = (n_a <= n_b) ? n_b : n_a; // n <-  max(length(a), length(b))
    SEXP r_ = PROTECT(allocVector(REALSXP, n));
    double *a = REAL(a_), *b = REAL(b_),
	*r = REAL(r_);
    for(i=0; i < n; i++) {
	r[i] = algdiv(a[i % n_a], b[i % n_b]);
    }
    UNPROTECT(1);
    return r_;
}
