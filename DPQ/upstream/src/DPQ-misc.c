/*
 *  Miscellanous  C  functions for pkg  'DPQ'
 *  ------
 */

// #include <float.h> /* DBL_MIN etc */

// for error() or if(verbose) ..
#include <R_ext/Print.h>

#include "DPQpkg.h"


// .Call()ed  from R --> ../R/utils.R <<<<<<<<

// Computes 'log(1 + X) - X'  accurately even for  |x| << 1
SEXP R_log1pmx(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n);
    const double *x = REAL(x_);
    double *r = REAL(r_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = log1pmx(x[i]); // R API log1pmx() [in R's src/nmath/pgamma.c]
    UNPROTECT(1);
    return r_;
}


// Computes 'log(1 + exp(X))' accurately, notably for large x, e.g., x > 720.
SEXP R_log1pexp(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n);
    const double *x = REAL(x_);
    double *r = REAL(r_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = log1pexp(x[i]);
    UNPROTECT(1);
    return r_;
}

#include <Rversion.h>
#if R_VERSION < R_Version(4,1,0)
// From --- src/nmath/{dpq.h, plogis.c} ---- :
// log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x)) :
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))
double log1mexp(double x) { return R_Log1_Exp(-x); }
#endif


/** Computes 'log(1 - exp(-x))' accurately, carefully for two regions of x,
 * optimally cutting off at log 2 (= 0.693147..),
 * if(x <= log(2)  log(-expm1(-x))  else  log1p(-exp(-x))
 */
SEXP R_log1mexp(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n);
    const double *x = REAL(x_);
    double *r = REAL(r_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = log1mexp(x[i]);
    UNPROTECT(1);
    return r_;
}


/* Computes 'log(gamma(X + 1))'  accurately even for small x, i.e., 0 < x < 0.5.
*  called from  DPQ :: lgamma1pC() */
SEXP R_lgamma1p(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP)));
    SEXP r_ = PROTECT(allocVector(REALSXP, n));
    const double *x = REAL(x_);
    double *r = REAL(r_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = lgamma1p(x[i]);
    UNPROTECT(2);
    return r_;
}

/* in C have     r = frexp (x, &e);
 * here, as "proper function" of _1_ argument, returning "two results" as list */
SEXP R_frexp(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP)));
    static const char *ans_nms[] = {"r", "e", ""};
    SEXP ans = PROTECT(mkNamed(VECSXP, ans_nms)), r_, e_; // --> ans = list(r = r_, e = e_)
    SET_VECTOR_ELT(ans, 0, r_ = PROTECT(allocVector(REALSXP, n)));
    SET_VECTOR_ELT(ans, 1, e_ = PROTECT(allocVector(INTSXP,  n)));

    const double *x = REAL(x_);
    double *r = REAL(r_);
    int *e = INTEGER(e_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = frexp(x[i], e+i);
    UNPROTECT(4);
    return ans;
}

// ldexp(f, E) := f * 2^E
SEXP R_ldexp(SEXP f_, SEXP E_)
{
    R_xlen_t
	n_f = XLENGTH(PROTECT(f_ = isReal   (f_) ? f_ : coerceVector(f_, REALSXP))),
	n_E = XLENGTH(PROTECT(E_ = isInteger(E_) ? E_ : coerceVector(E_, INTSXP))),
	n = (!n_f || !n_E) ? 0 : (n_f >= n_E) ? n_f : n_E; // and recycle f_ or E_ when needed
    SEXP r_ = PROTECT(allocVector(REALSXP, n));
    if(n) {
	const double *f = REAL(f_);
	double *r = REAL(r_);
	int *E = INTEGER(E_);
	for(R_xlen_t i=0; i < n; i++)
	    r[i] = ldexp(f[i % n_f], E[i % n_E]);
    }
    UNPROTECT(3);
    return r_;
}

/* in C have     fr = modf (x, &i);
 * modf stores the integer part in *i, and returns the fractional part.
 *
 * returning the "two results" as list(fr=., i=.) in R */
SEXP R_modf(SEXP x_)
{
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP)));
    static const char *ans_nms[] = {"fr", "i", ""};
    SEXP ans = PROTECT(mkNamed(VECSXP, ans_nms)), r_, i_; // --> ans = list(fr = r_, i = i_)
    SET_VECTOR_ELT(ans, 0, r_ = PROTECT(allocVector(REALSXP, n)));
    SET_VECTOR_ELT(ans, 1, i_ = PROTECT(allocVector(REALSXP, n)));
    const double *x = REAL(x_);
    double *r = REAL(r_), *e = REAL(i_);
    for(R_xlen_t i=0; i < n; i++)
	r[i] = modf(x[i], e+i);
    UNPROTECT(4);
    return ans;
}

/* Compute x_ ^ y_  */
SEXP dpq_pow(SEXP x_, SEXP y_)
{
    R_xlen_t
	nx = XLENGTH(PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP))),
	ny = XLENGTH(PROTECT(y_ = isReal(y_) ? y_ : coerceVector(y_, REALSXP))),
	n = (!nx || !ny) ? 0 : (nx <= ny) ? ny : nx;
    SEXP r_ = PROTECT(allocVector(REALSXP, n));
    if(n) {
	const double *x = REAL(x_), *y = REAL(y_);
	double *r = REAL(r_);
	for(R_xlen_t i=0; i < n; i++)
	    r[i] = R_pow(x[i % nx], y[i % ny]);
    }
    UNPROTECT(3);
    return r_;
}

/* Compute x_ ^ y_ -- y integer */
SEXP dpq_pow_di(SEXP x_, SEXP y_)
{
    R_xlen_t
	nx = XLENGTH(PROTECT(x_ = isReal   (x_) ? x_ : coerceVector(x_, REALSXP))),
	ny = XLENGTH(PROTECT(y_ = isInteger(y_) ? y_ : coerceVector(y_,  INTSXP))),
	n = (!nx || !ny) ? 0 : (nx <= ny) ? ny : nx;
    SEXP r_ = PROTECT(allocVector(REALSXP, n));
    if(n) {
	const double *x = REAL(x_);
	double *r = REAL(r_);
	const int *y = INTEGER(y_);
	for(R_xlen_t i=0; i < n; i++)
	    r[i] = R_pow_di(x[i % nx], y[i % ny]);
    }
    UNPROTECT(3);
    return r_;
}

