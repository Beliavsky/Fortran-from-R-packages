#include "DPQpkg.h"

#undef min
#define min(a,b) ((a < b)?a:b)
#undef max
#define max(a,b) ((a > b)?a:b)

/*---------------------------------------------vvvvvvvvvvvvvvvvvvvvv-------
 * Cut'n'paste from R's sources  ~/R/D/r-devel/R/src/nmath/toms708.c
 * [as of 2022-03-19, lines  512-659]
 * replacing s/betaln/lbeta/ in one place
 *-------------------------------------------------------------------------*/

double bpser(double a, double b, double x, double eps, int *err_bp, int log_p, Rboolean verbose)
{
/* -----------------------------------------------------------------------
 * Power SERies expansion for evaluating I_x(a,b) when
 *	       b <= 1 or b*x <= 0.7.   eps is the tolerance used.
 * NB: if log_p is TRUE, also use it if   (b < 40  & lambda > 650)   where
 *                       lambda := a y - b x = (a + b)y - b  = a - (a+b)x   {x + y == 1}
 * ----------------------------------------------------------------------- */

    if (x == 0.) {
	return R_D__0;
    }

    /* for DPQ: */
    if (x == 1. || a == 0.) return R_D__1;

/* ----------------------------------------------------------------------- */
/*	      compute the factor  x^a/(a*Beta(a,b)) */
/* ----------------------------------------------------------------------- */
    double ans, c, z, a0 = min(a,b);

    if (a0 >= 1.) { /*		------	 1 <= a0 <= b0  ------ */
	z = a * log(x) - lbeta(a, b);
	ans = log_p ? z - log(a) : exp(z) / a;
    }
    else { /*  0 <= a0 < 1  */
	double t, u, apb, b0 = max(a,b);
	if (b0 < 8.) {

	    if (b0 <= 1.) { /*	------  a0 < 1  and  a0 <= b0 <= 1  ------ */

		if(log_p) {
		    ans = a * log(x);
		} else {
		    ans = pow(x, a);
		    if (ans == 0.) { /* once underflow, always underflow .. */
			if(verbose) REprintf(" bpser(a=%g, b=%g, x=%g): x^a underflows to 0",
					 a,b,x);
			return ans;
		    }
		}
		apb = a + b;
		if (apb > 1.) {
		    u = a + b - 1.;
		    z = (gam1(u, TRUE, verbose >= 2) + 1.) / apb;
		} else {
		    z = gam1(apb, TRUE, verbose >= 2) + 1.;
		}
		c = (gam1(a, TRUE, verbose >= 2) + 1.) * (gam1(b, TRUE, verbose >= 2) + 1.) / z;

		if(log_p) /* FIXME ? -- improve quite a bit for c ~= 1 */
		    ans += log(c * (b / apb));
		else
		    ans *=  c * (b / apb);

	    } else { /* 	------	a0 < 1 < b0 < 8	 ------ */

		u = gamln1(a0, TRUE);
		int m = (int)(b0 - 1.);
		if (m >= 1) {
		    c = 1.;
		    for (int i = 1; i <= m; ++i) {
			b0 += -1.;
			c *= b0 / (a0 + b0);
		    }
		    u += log(c);
		}

		z = a * log(x) - u;
		b0 += -1.; // => b0 in (0, 7)
		apb = a0 + b0;
		if (apb > 1.) {
		    u = a0 + b0 - 1.;
		    t = (gam1(u, TRUE, verbose >= 2) + 1.) / apb;
		} else {
		    t = gam1(apb, TRUE, verbose >= 2) + 1.;
		}

		if(log_p) /* FIXME? potential for improving log(t) */
		    ans = z + log(a0 / a) + log1p(gam1(b0, TRUE, verbose >= 2)) - log(t);
		else
		    ans = exp(z) * (a0 / a) * (gam1(b0, TRUE, verbose >= 2) + 1.) / t;
	    }

	} else { /* 		------  a0 < 1 < 8 <= b0  ------ */

	    u = gamln1(a0, TRUE) + algdiv(a0, b0);
	    z = a * log(x) - u;

	    if(log_p)
		ans = z + log(a0 / a);
	    else
		ans = a0 / a * exp(z);
	}
    }
    if(verbose) REprintf(" bpser(a=%g, b=%g, x=%g, log=%d, eps=%g): %s = %.14g;",
		     a,b,x, log_p, eps,
		     log_p ? "log(x^a/(a*B(a,b)))" : "x^a/(a*B(a,b))",  ans);
    if (ans == R_D__0 || (!log_p && a <= eps * 0.1)) {
	if(verbose) REprintf(" = final answer\n");
	return ans;
    }
    else if(verbose) REprintf("\n");

/* ----------------------------------------------------------------------- */
/*		       COMPUTE THE SERIES */
/* ----------------------------------------------------------------------- */
    double tol = eps / a,
	n = 0.,
	sum = 0., w;
    c = 1.;
    do { // sum is alternating as long as n < b (<==> 1 - b/n < 0)
	n += 1.;
	c *= (0.5 - b / n + 0.5) * x;
	w = c / (a + n);
	sum += w;
    } while (n < 1e7 && fabs(w) > tol);
    if(fabs(w) > tol) { // the series did not converge (in time)
	// warn only when the result seems to matter:
	if(( log_p && !(a*sum > -1. && fabs(log1p(a * sum)) < eps*fabs(ans))) ||
	   (!log_p && fabs(a*sum + 1.) != 1.)) {
	    if(*err_bp >= 0) // caller can specify err_bp = -1  to suppress this warning
	    MATHLIB_WARNING5(
		" bpser(a=%g, b=%g, x=%g,...) did not converge (n=1e7, |w|/tol=%g > 1; A=%g)",
		a,b,x, fabs(w)/tol, ans);
	    *err_bp = 1;
	}
    }
    if(verbose) REprintf("  -> n=%.0f iterations, |w|=%g %s %g=tol:=eps/a ==> a*sum=%g %s -1\n",
		     n, fabs(w), (fabs(w) > tol) ? ">!!>" : "<=", tol,
		     a*sum, (a*sum > -1.) ? ">" : "<=");
    if(log_p) {
	if (a*sum > -1.) ans += log1p(a * sum);
	else {
	    if(ans > ML_NEGINF) {
		if(*err_bp >= -1) // caller can specify err_bp = -2  to suppress both warnings
		MATHLIB_WARNING3(
		    "pbeta(*, log.p=TRUE) -> bpser(a=%g, b=%g, x=%g,...) underflow to -Inf",
		    a,b,x);
		*err_bp = 2;
	    }
	    // FIXME ? rather keep first order term ans = log(x^a/(a*B(a,b))) from above
	    ans = ML_NEGINF;
	}
    } else if (a*sum > -1.)
	ans *= (a * sum + 1.);
    else // underflow to
	ans = 0.;
    return ans;
} /* bpser */
// ---------  ------------------------- end cut'n'paste from R sources

// To be  .Call()ed  from R :
SEXP R_bpser(SEXP a_, SEXP b_, SEXP x_, SEXP eps_, SEXP log_p_, SEXP verbose_, SEXP warn_)
{
    PROTECT(a_ = isReal(a_) ? a_ : coerceVector(a_, REALSXP));
    PROTECT(b_ = isReal(b_) ? b_ : coerceVector(b_, REALSXP));
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    R_xlen_t i,
	n_a = XLENGTH(a_),
	n_b = XLENGTH(b_),
	n_x = XLENGTH(x_), n;
    if(!n_x || !n_a || !n_b)
	n = 0;
    else { // n <-  max(length(a), length(b), length(x))
	n = (n_a <= n_b) ? n_b : n_a;
	if(n_x > n) n = n_x;
    }
    static const char *ans_nms[] = {"r", "err", ""};
    SEXP ans = PROTECT(mkNamed(VECSXP, ans_nms)), r_, ier_; // --> ans = list(r = r_, err = ier_)
    SET_VECTOR_ELT(ans, 0, r_   = PROTECT(allocVector(REALSXP, n)));
    SET_VECTOR_ELT(ans, 1, ier_ = PROTECT(allocVector(INTSXP,  n)));
    if(n >= 1) {
	int *ierr = INTEGER(ier_);
	double *a = REAL(a_), *b = REAL(b_), *x = REAL(x_), eps = asReal(eps_),
	    *r = REAL(r_);
	Rboolean verbose = asLogical(verbose_),
	    warn = asLogical(warn_);
	int log_p = (int)asLogical(log_p_);

	for(i=0; i < n; i++) {
	    ierr[i] = warn ? 0 : -2;
	    r[i] = bpser(a[i % n_a], b[i % n_b], x[i % n_x], eps, &ierr[i], log_p, verbose);
	}
    }
    UNPROTECT(6);
    return ans;
}

#ifdef DEBUG_gam1
# define R_ifDEBUG_printf(...) REprintf(__VA_ARGS__)
#else
# define R_ifDEBUG_printf(...)
#endif



// == 1/gamma(a+1) - 1   {not clear why this is needed}
double gam1(double a, Rboolean warn_if, Rboolean verbose)
{
/*     ------------------------------------------------------------------ */
/*     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 <= A <= 1.5 */
/*     ------------------------------------------------------------------ */

    double t = a,
	d = a - 0.5;
    // t := if(a > 1/2)  a-1  else  a
    if (d > 0.)
	t = d - 0.5;

    double w, bot, top;
    if (t < 0.) { /* L30: */
	static double
	    r[9] = { -.422784335098468,-.771330383816272,
		     -.244757765222226,.118378989872749,9.30357293360349e-4,
		     -.0118290993445146,.00223047661158249,2.66505979058923e-4,
		     -1.32674909766242e-4 },
	    s1 = .273076135303957,
	    s2 = .0559398236957378;

	top = (((((((r[8] * t + r[7]) * t + r[6]) * t + r[5]) * t + r[4]
		     ) * t + r[3]) * t + r[2]) * t + r[1]) * t + r[0];
	bot = (s2 * t + s1) * t + 1.;
	w = top / bot;
	if(verbose) REprintf("  gam1(a = %.15g): t < 0: w=%.15g\n", a, w);
	if(warn_if && a > 1.5)
	    warning("gam1(a = %g) is used outside of documented boundaries, [-0.5, 1.5]\n",
		    a);

	if (d > 0.)
	    return t * w / a;
	else
	    return a * (w + 0.5 + 0.5);

    } else if (t == 0) { // L10: a in {0, 1}
	return 0.;

    } else { /* t > 0;  L20: */
	static double
	    p[7] = { .577215664901533,-.409078193005776,
		     -.230975380857675,.0597275330452234,.0076696818164949,
		     -.00514889771323592,5.89597428611429e-4 },
	    q[5] = { 1.,.427569613095214,.158451672430138,
		     .0261132021441447,.00423244297896961 };

	top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]
		   ) * t + p[1]) * t + p[0];
	bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.;
	w = top / bot;

	if(verbose) REprintf("  gam1(a = %.15g): t > 0: (is a < 1.5 ?)  w=%.15g\n",
			     a, w);
	if(warn_if && a < -0.5)
	    warning("gam1(a = %g) is used outside of documented boundaries, [-0.5, 1.5]\n",
		    a);
	if (d > 0.) /* L21: */
	    return t / a * (w - 0.5 - 0.5);
	else
	    return a * w;
    }
} /* gam1 */

// == lgamma(a+1) --- bpser() used it only for 0 <= a0 < 1
double gamln1(double a, Rboolean warn_if)
{
/* ----------------------------------------------------------------------- */
/*     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25 */
/* ----------------------------------------------------------------------- */

    double w;
    if (a < 0.6) {
	if(warn_if && a < -0.2)
	    warning("gamln1(a = %g) is used outside of documented boundaries, [-0.2, 1.25]\n",
		    a);
	static double p0 = .577215664901533;
	static double p1 = .844203922187225;
	static double p2 = -.168860593646662;
	static double p3 = -.780427615533591;
	static double p4 = -.402055799310489;
	static double p5 = -.0673562214325671;
	static double p6 = -.00271935708322958;
	static double q1 = 2.88743195473681;
	static double q2 = 3.12755088914843;
	static double q3 = 1.56875193295039;
	static double q4 = .361951990101499;
	static double q5 = .0325038868253937;
	static double q6 = 6.67465618796164e-4;
	w = ((((((p6 * a + p5)* a + p4)* a + p3)* a + p2)* a + p1)* a + p0) /
	    ((((((q6 * a + q5)* a + q4)* a + q3)* a + q2)* a + q1)* a + 1.);
	return -(a) * w;
    }
    else { /* 0.6 <= a <= 1.25 */
	if(warn_if && a > 1.25)
	    warning("gamln1(a = %g) is used outside of documented boundaries, [-0.2, 1.25]\n",
		    a);
	static double r0 = .422784335098467;
	static double r1 = .848044614534529;
	static double r2 = .565221050691933;
	static double r3 = .156513060486551;
	static double r4 = .017050248402265;
	static double r5 = 4.97958207639485e-4;
	static double s1 = 1.24313399877507;
	static double s2 = .548042109832463;
	static double s3 = .10155218743983;
	static double s4 = .00713309612391;
	static double s5 = 1.16165475989616e-4;
	double x = a - 0.5 - 0.5;
	w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) /
	    (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x + 1.);
	return x * w;
    }
} /* gamln1 */

//------ gam1()  &  gamln1() ---- R interface --- to assess their accuracy, etc

SEXP R_gam1 (SEXP a_, SEXP warn_if_, SEXP verbose_) {
    PROTECT(a_ = isReal(a_) ? a_ : coerceVector(a_, REALSXP));
    R_xlen_t i, n = XLENGTH(a_);
    Rboolean
	warn_if = asLogical(warn_if_),
	verbose = asLogical(verbose_);
    SEXP res = PROTECT(allocVector(REALSXP, n));
    double *a = REAL(a_), *r = REAL(res);
    for(i=0; i < n; i++) {
	r[i] = gam1(a[i], warn_if, verbose);
    }
    UNPROTECT(2);
    return res;
}

SEXP R_gamln1 (SEXP a_, SEXP warn_if_) {
    PROTECT(a_ = isReal(a_) ? a_ : coerceVector(a_, REALSXP));
    Rboolean warn_if = asLogical(warn_if_);
    R_xlen_t i, n = XLENGTH(a_);
    SEXP res = PROTECT(allocVector(REALSXP, n));
    double *a = REAL(a_), *r = REAL(res);
    for(i=0; i < n; i++) {
	r[i] = gamln1(a[i], warn_if);
    }
    UNPROTECT(2);
    return res;
}

