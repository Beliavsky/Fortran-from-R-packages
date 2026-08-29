#include "DPQpkg.h"

/* calling dpsifn() from R's  <Rsource>/src/nmath/polygamma.c :

 * From R, currently only used for kode = 1, m = 1 :

 void dpsifn(double x, int n, int kode, int m, double *ans, int *nz, int *ierr)
*/

SEXP R_dpsifn(SEXP x_, SEXP m_, SEXP deriv_1, SEXP kode2_) {
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    // other args are checked via  asReal(.), asInteger() below
    R_xlen_t i;
    if(XLENGTH(x_) > INT_MAX)
	error("length(%s) = %lld > max.int = %d", "x", (long long)XLENGTH(x_), INT_MAX);
    int n  = (int)XLENGTH(x_),
	m  = asInteger(m_),
	d1 = asInteger(deriv_1);
    /* result: a list(dpsi = <m x n> matrix of dpsi() values,
     *                nz   =  n-vector,
     *                ierr =  n-vector)
     */
    const char *nms[] = {"dpsi", "nz", "ierr", ""};
    SEXP ans = PROTECT(Rf_mkNamed(VECSXP, nms));
    SET_VECTOR_ELT(ans, 0, allocMatrix(REALSXP, m, n));// 'dpsi'
    SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, n));   // 'nz'
    SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, n));   // 'ierr'

    double *x = REAL(x_),
	*r = REAL(VECTOR_ELT(ans, 0));
    int *nz   = INTEGER(VECTOR_ELT(ans, 1)),
	*ierr = INTEGER(VECTOR_ELT(ans, 2));
    Rboolean k2 = asLogical(kode2_);
    //  R[_i_,_j_] := r[m*(_i_) + _j_]
    for(i=0; i < n; i++) {
	dpsifn(x[i], /* n = */ d1, /* kode = */ k2 ? 2 : 1, m,
	       /* *ans = */ r+m*i,
	       nz+i,
	       ierr+i);
    }
    UNPROTECT(2);
    return ans;
}
