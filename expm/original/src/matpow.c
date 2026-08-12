/*
 *  Power of a matrix x^k := x x ... x, where x in an (n x n) matrix
 *  and k is an *integer*. Based on code originally written by Vincent
 *  Goulet for package actuar (and inspired from Octave in file
 *  .../src/xpow.cc) with slight shortcuts by Martin Maechler:
 */

#include "matpow.h"

/* .Call() this from R : */
SEXP R_matpow(SEXP x, SEXP k)
{
    if(!isMatrix(x)) {
	error(_("not a matrix"));
	/*-Wall */ return R_NilValue;
    }
    else {
	SEXP dims = getAttrib(x, R_DimSymbol);
	int n = INTEGER(dims)[0],
	    k_ = INTEGER(k)[0]; /* need copy, as it is altered in matpow() */

	if (n != INTEGER(dims)[1])
	    error(_("non-square matrix"));
	if (n == 0)
	    return(allocMatrix(REALSXP, 0, 0));

	SEXP z, x_ = duplicate(x); PROTECT_INDEX xpi;
	PROTECT_WITH_INDEX(x_, &xpi);
	if(isComplex(x)) {
	    z = PROTECT(allocMatrix(CPLXSXP, n, n));
	    matpow_z(COMPLEX(x_), n, k_, COMPLEX(z));
	} else { // real / double
	    if (!isReal(x)) /* coercion to numeric */
		REPROTECT(x_ = coerceVector(x_, REALSXP), xpi);
	    z = PROTECT(allocMatrix(REALSXP, n, n));
	    matpow(REAL(x_), n, k_, REAL(z));
	}
	setAttrib(z, R_DimNamesSymbol,
		  getAttrib(x, R_DimNamesSymbol));
	UNPROTECT(2);
	return z;
    }
}


/* Compute z := x %^% k, x an (n x n) square "matrix" in column-order;
 * NB: x will be altered! The caller must make a copy if needed */
void matpow(double *x, int n, int k, double *z)
{
    if (k == 0) { /* return identity matrix */
	int i, j;
	for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++)
		z[i * n + j] = (i == j) ? 1.0 : 0.0;
	return;
    }
    else if (k < 0) {
	error(_("power must be a positive integer; use solve() directly for negative powers"));
    }
    else { /* k >= 1 */
	static const char *transa = "N";
	static const double one = 1.0, zero = 0.0;
	size_t nSqr = n * ((size_t) n);
	double /* temporary matrix */
	    *tmp  = (double *) R_alloc(nSqr, sizeof(double));

	/* Take powers in multiples of 2 until there is only one
	 * product left to make. That is, if k = 5, compute (x * x),
	 * then ((x * x) * (x * x)) and finally ((x * x) * (x * x)) * x.
	 */
	Memcpy(z, x, nSqr);

	k--;
	while (k > 0) {
	    if (k & 1) {	/* z := z * x */
		F77_CALL(dgemm)(transa, transa, &n, &n, &n, &one,
				z, &n, x, &n, &zero, tmp, &n FCONE FCONE);
		Memcpy(z, tmp, nSqr);
	    }
	    if(k == 1)
		break;
	    k >>= 1; /* efficient division by 2; now have k >= 1 */

	    /* x := x * x */
	    F77_CALL(dgemm)(transa, transa, &n, &n, &n, &one,
			    x, &n, x, &n, &zero, tmp, &n FCONE FCONE);
	    Memcpy(x, tmp, nSqr);
	}
    }
}

/* Compute z := x %^% k, x an (n x n) square "matrix" in column-order;
 * NB: x will be altered! The caller must make a copy if needed */
void matpow_z(Rcomplex *x, int n, int k, Rcomplex *z)
{
    if (k == 0) { /* return identity matrix */
	for (int i = 0; i < n; i++)
	    for (int j = 0; j < n; j++) {
		z[i * n + j].r = (i == j) ? 1. : 0.;
		z[i * n + j].i = 0.;
	    }
	return;
    }
    else if (k < 0) {
	error(_("power must be a positive integer; use solve() directly for negative powers"));
    }
    else { /* k >= 1 */
	static const char *transa = "N";
	static const Rcomplex
	    zone  = { .r = 1., .i = 0.},
	    zzero = { .r = 0., .i = 0.};
	size_t nSqr = n * ((size_t) n);
	Rcomplex /* temporary matrix */
	    *tmp  = (Rcomplex *) R_alloc(nSqr, sizeof(Rcomplex));

	/* Take powers in multiples of 2 until there is only one
	 * product left to make. That is, if k = 5, compute (x * x),
	 * then ((x * x) * (x * x)) and finally ((x * x) * (x * x)) * x.
	 */
	Memcpy(z, x, nSqr);

	k--;
	while (k > 0) {
	    if (k & 1) {	/* z := z * x */
		F77_CALL(zgemm)(transa, transa, &n, &n, &n, &zone,
				z, &n, x, &n, &zzero, tmp, &n FCONE FCONE);
		Memcpy(z, tmp, nSqr);
	    }
	    if(k == 1)
		break;
	    k >>= 1; /* efficient division by 2; now have k >= 1 */

	    /* x := x * x */
	    F77_CALL(zgemm)(transa, transa, &n, &n, &n, &zone,
			    x, &n, x, &n, &zzero, tmp, &n FCONE FCONE);
	    Memcpy(x, tmp, nSqr);
	}
    }
}
