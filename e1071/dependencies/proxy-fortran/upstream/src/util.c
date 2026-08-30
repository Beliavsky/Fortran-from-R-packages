#include <R.h>
#include <Rdefines.h>

// arrayIndex.c
extern SEXP _int_array_subscript(int, SEXP, const char *, const char *, SEXP, Rboolean, SEXP);

// subset a dist object. in order to preserve symmetry
// we allow only one subset index.
//
// notes: (1) by definition we return a zero-length vector
// for subscripts of length less than two. (2) coercing to
// real is convenient but slightly inefficient.
//
// ceeboo 2006, 2007

SEXP R_subset_dist(SEXP R_x, SEXP s) {
    if (!inherits(R_x, "dist"))
	error("'x' not of class dist");
    int i, j, k, si, sj, nx, ns;
    SEXP x = R_x, r, d;
    R_xlen_t l;

    nx = 1 + (int) sqrt(2*XLENGTH(x));
    l = nx;
    if (XLENGTH(x) != l*(l-1)/2)
	error("'x' invalid length");
    
    if (TYPEOF(x) != REALSXP) 
	PROTECT(x = coerceVector(R_x, REALSXP));

    PROTECT(r = allocArray(INTSXP, PROTECT(ScalarInteger(0))));
    UNPROTECT(1);
    INTEGER(getAttrib(r, R_DimSymbol))[0] = nx;

    d = getAttrib(x, install("Labels"));
    if (!isNull(d)) {
	SEXP t;

	if (TYPEOF(d) != STRSXP)
	    error("'Labels' not of type character");
	if (LENGTH(d) != nx)
	    error("'Labels' invalid length");
	setAttrib(r, R_DimNamesSymbol, PROTECT(t = allocVector(VECSXP, 1)));
	UNPROTECT(1);
	SET_VECTOR_ELT(t, 0, d);
    }

#ifdef _COMPAT_
    PROTECT(s = arraySubscript(0, s, GET_DIM(r), getAttrib, (STRING_ELT), r));
#else
    PROTECT(s = _int_array_subscript(0, s, "dim", "dimnames", r,
						  TRUE, R_NilValue));
#endif
    ns = LENGTH(s);
    
    for (k = 0; k < ns; k++)
	if (INTEGER(s)[k] == NA_INTEGER)
	    error("'s' invalid subscript(s)");
	else
	    INTEGER(s)[k]--;

    l = ns;
    PROTECT(r = allocVector(REALSXP, l*(l-1)/2));

    l = 0;
    for (i = 0; i < ns-1; i++) {
	si = INTEGER(s)[i];
	for (j = i+1; j < ns; j++) {
	    sj = INTEGER(s)[j];
	    if (si == sj)
		REAL(r)[l++] = NA_REAL;
	    else 
		REAL(r)[l++] = 
		(si > sj)    ? REAL(x)[si+sj*(nx-1)-sj*(sj+1)/2-1]
			     : REAL(x)[sj+si*(nx-1)-si*(si+1)/2-1];
	}
	R_CheckUserInterrupt();
    }

    if (x == R_x)
	copyMostAttrib(R_x, r);

    setAttrib(r, install("Size"), PROTECT(ScalarInteger(ns)));
    UNPROTECT(1);
    if (!isNull(d)) {
	SEXP t;
	
	setAttrib(r, install("Labels"), PROTECT(t = allocVector(STRSXP, ns)));
	UNPROTECT(1);
	for (k = 0; k < ns; k++)
	    SET_STRING_ELT(t, k, STRING_ELT(d, INTEGER(s)[k]));
    }
    
    UNPROTECT(3);
    if (x != R_x)
	UNPROTECT(1);
    
    return r;
}

// compute the rowSums for an R dist object. due to
// symmetry this is equivalent to colSums. rowMeans
// are not implemented as these can be easily obtained
// from the values of rowSums.
//
// na_rm implements the usual meaning of omitting NA
// and NaN values, where the sum of the empty set is
// zero.
//
// ceeboo 2006, 2007

SEXP R_rowSums_dist(SEXP R_x, SEXP na_rm) {
    if (!inherits(R_x, "dist"))
	error("'x' not of class dist");
    if (isNull(na_rm) || TYPEOF(na_rm) != LGLSXP)
        error("'na.rm' not of type logical");
    int i, j, n;
    SEXP x = R_x, r;
    R_xlen_t k;

    n = 1 + (int) sqrt(2*XLENGTH(x));
    k = n;
    if (XLENGTH(x) != k*(k-1)/2)
        error("'x' invalid length");
   
    if (TYPEOF(x) != REALSXP) 
	PROTECT(x = coerceVector(R_x, REALSXP));

    PROTECT(r = allocVector(REALSXP, n));

    memset(REAL(r), 0, sizeof(double)*n);

    k = 0;
    for (i = 0; i < n-1; i++) {
	for (j = i+1; j < n; j++) {
	    double z = REAL(x)[k++];
	    if (!R_FINITE(z)) {
		if (ISNAN(z)) {
		    if (LOGICAL(na_rm)[0] == TRUE)
			continue;
		    REAL(r)[i] = REAL(r)[j] = (ISNA(z)) ? NA_REAL : R_NaN;
		} else
		    REAL(r)[i] = REAL(r)[j] = z;
		break;
            }
            REAL(r)[i] += z;
            REAL(r)[j] += z;
        }
	R_CheckUserInterrupt();
    }
    setAttrib(r, R_NamesSymbol, getAttrib(x, install("Labels")));
    
    UNPROTECT(1);
    if (x != R_x)
	UNPROTECT(1);

    return r;
}

// produce row or column indexes

SEXP R_row_dist(SEXP x, SEXP col) {
    if (!inherits(x, "dist"))
	error("'x' not of class dist");
    if (isNull(col) || TYPEOF(col) != LGLSXP)
	error("'col' not of type logical");
    int i, j, nx;
    SEXP r;
    R_xlen_t l;

    nx = 1 + (int) sqrt(2*XLENGTH(x));
    l = nx;
    if (XLENGTH(x) != l*(l-1)/2)
	error("'x' invalid length");

    PROTECT(r = allocVector(INTSXP, XLENGTH(x)));
    
    l = 0;
    for (j = 1; j < nx; j++)
	for (i = j+1; i < nx+1; i++) 
	    INTEGER(r)[l++] = (*LOGICAL(col)) ? j : i;

    UNPROTECT(1);

    return r;
}

//
