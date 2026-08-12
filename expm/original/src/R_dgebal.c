#include <ctype.h>
	/* strlen(), toupper() .. */

#include "expm.h"

static char ebal_type(const char *typstr)
{
    if (strlen(typstr) != 1)
	error(_("argument type='%s' must be a character string of string length 1"),
	      typstr);
    char typup = (char)toupper((int)*typstr);
    /* if (typup == '1') typup = 'O'; /\* alias *\/ */
    if (typup != 'N' && typup != 'P' && typup != 'S' && typup != 'B')
	error(_("argument type='%s' must be one of 'N', 'P', 'S', or 'B'"),
	      typstr);
    return typup;
}

SEXP R_dgebal(SEXP x, SEXP type)
{
    SEXP dims, z, Scale, i_1, i_2, ans, nms;
    char typnm[] = {'\0', '\0'}; // only the first is changed; 2nd = final \0 string terminator
    int n, info, nprot = 2;

    if (!isNumeric(x) || !isMatrix(x)) // isNum : integer, logical, or double ("real")
	error(_("invalid 'x': not a numeric (classical R) matrix"));
    dims = getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    if (n != INTEGER(dims)[1])
	error(_("non-square matrix"));
    typnm[0] = ebal_type(CHAR(asChar(type)));
    if (isInteger(x) || isLogical(x)) {
	nprot++;
	x = PROTECT(coerceVector(x, REALSXP));
    }
    else if(n > 0 && typnm[0] == 'S') {
	/* FIXME: if 'x' contains +/- Inf dgebal() loops infinitely <==> LAPACK "bug"
	   ----- fix in ..../R/src/modules/lapack/dlapack.f.~dgebal-Inf-patch~
	   But that does *not* help for external Lapack libraries */
	double *dx = REAL(x), aMax = 0.; // aMax := max_{i,j} |x[i,j]|
	for(int i=0; i < n*n; i++)
	    if(aMax < dx[i]) aMax = dx[i];
	if(aMax == R_PosInf)
	    error(_("R_dgebal(*, type=\"S\"): Infinite matrix entry"));
    }

    PROTECT(ans = allocVector(VECSXP, 4));
    PROTECT(nms = allocVector(STRSXP, 4));

    SET_STRING_ELT(nms, 0, mkChar("z"));
    SET_VECTOR_ELT(ans, 0, (z = duplicate(x)));

    /* permutation or scale array */
    SET_STRING_ELT(nms, 1, mkChar("scale"));
    SET_VECTOR_ELT(ans, 1, (Scale = allocVector(REALSXP, n)));

    SET_STRING_ELT(nms, 2, mkChar("i1"));
    SET_VECTOR_ELT(ans, 2, (i_1 = allocVector(INTSXP, 1)));

    SET_STRING_ELT(nms, 3, mkChar("i2"));
    SET_VECTOR_ELT(ans, 3, (i_2 = allocVector(INTSXP, 1)));

    if(n > 0) {
	F77_CALL(dgebal)(typnm, &n, REAL(z), &n, INTEGER(i_1), INTEGER(i_2),
			 REAL(Scale), &info FCONE);
	if (info)
	    error(_("LAPACK's dgebal(%s) returned info code %d"), typnm, info);
    }

    setAttrib(ans, R_NamesSymbol, nms);
    /* now return  list(z, scale[], i1, i2) */
    UNPROTECT(nprot);
    return ans;
} // R_dgebal()

//-------------- Now for complex matrices ----------------------------

#include <complex.h>
#include <R_ext/Complex.h>

//  CMPLX(.) "must be" in <complex.h> with 'C11' and newer, still -- fails on Winbuilder (2024-08-09)
/* # ifndef CMPLX */
/* # define CMPLX(x, y)  ((double complex)((double)(x) + _Imaginary_I * (double)(y))) */
/* # endif */
// rather for now (need only cabs(.)), this seems safer :
#ifdef CMPLX
static R_INLINE double R_cabs(Rcomplex z) { return cabs(CMPLX(z.r, z.i)); }
#else // currently necessary on Winbuilder
static R_INLINE double R_cabs(Rcomplex z) { return hypot(z.r, z.i); }
#endif

SEXP R_zgebal(SEXP x, SEXP type)
{
    SEXP dims, z, Scale, i_1, i_2, ans, nms;
    char typnm[] = {'\0', '\0'}; // only the first is changed; 2nd = final \0 string terminator
    int n, info, nprot = 2;

    if (!isComplex(x) || !isMatrix(x))
	error(_("invalid 'x': not a complex (classical R) matrix"));
    dims = getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    if (n != INTEGER(dims)[1])
	error(_("non-square matrix"));
    typnm[0] = ebal_type(CHAR(asChar(type)));
    if(n > 0 && typnm[0] == 'S') {
	/* FIXME: if 'x' contains +/- Inf dgebal() loops infinitely <==> LAPACK "bug"
	   ----- fix in ..../R/src/modules/lapack/dlapack.f.~dgebal-Inf-patch~
	   But that does *not* help for external Lapack libraries */
	Rcomplex *dx = COMPLEX(x);
	double aMax = 0.; // aMax := max_{i,j} |x[i,j]|
	for(int i=0; i < n*n; i++) {
	    double ai = R_cabs(dx[i]);
	    if(aMax < ai) aMax = ai;
	}
	if(aMax == R_PosInf)
	    error(_("R_zgebal(*, type=\"S\"): Infinite matrix entry"));
    }

    PROTECT(ans = allocVector(VECSXP, 4));
    PROTECT(nms = allocVector(STRSXP, 4));

    SET_STRING_ELT(nms, 0, mkChar("z"));
    SET_VECTOR_ELT(ans, 0, (z = duplicate(x)));

    /* permutation or scale array */
    SET_STRING_ELT(nms, 1, mkChar("scale"));
    SET_VECTOR_ELT(ans, 1, (Scale = allocVector(REALSXP, n)));

    SET_STRING_ELT(nms, 2, mkChar("i1"));
    SET_VECTOR_ELT(ans, 2, (i_1 = allocVector(INTSXP, 1)));

    SET_STRING_ELT(nms, 3, mkChar("i2"));
    SET_VECTOR_ELT(ans, 3, (i_2 = allocVector(INTSXP, 1)));

    if(n > 0) {
	F77_CALL(zgebal)(typnm, &n, COMPLEX(z), &n, INTEGER(i_1), INTEGER(i_2),
			 REAL(Scale), &info FCONE);
	if (info)
	    error(_("LAPACK's zgebal(%s) returned info code %d"), typnm, info);
    }

    setAttrib(ans, R_NamesSymbol, nms);
    /* now return  list(z, scale[], i1, i2) */
    UNPROTECT(nprot);
    return ans;
} // R_zgebal()
