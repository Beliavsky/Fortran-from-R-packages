
#include <Rdefines.h>

// wrapper functions for distance computation with
// user-supplied functions given common data types,
// such as matrix, list, and data.frames.
//
// note that the code is prone to breaking with new
// releases of R. therefore, always check against the
// reference examples.

/* compute auto- or cross-distances with a user-supplied
 * function given matrix data.
 * 
 * ceeboo 2006, 2007
 */

SEXP R_apply_dist_matrix(SEXP p) {
    int i, j, k, n, nx, ny, nz, m = 0;
    SEXP r, c, tx, ty;
    SEXP R_x, x, R_y, y, R_d, f;
    R_xlen_t l;

    p = CDR(p);
    if (length(p) < 4)
        error("invalid number of arguments");
    R_x = x = CAR(p); R_y = y = CADR(p);
    if (!isMatrix(x) || (!isNull(y) && !isMatrix(y)))
        error("invalid data parameter(s)");
    p = CDDR(p); R_d = CAR(p);
    if (TYPEOF(R_d) != LGLSXP)
	error("invalid option parameter");
    p = CDR(p); f = CAR(p); 
    if (!isFunction(f))
        error("invalid function parameter");
    p = CDR(p);

    if (isNull(y))
        y = x;  
    else
    if (LOGICAL(R_d)[0] == TRUE)
	m = 2;
    else
	m = 1;
    
    if ((n = INTEGER(GET_DIM(x))[1]) != INTEGER(GET_DIM(y))[1])
        error("the number of columns of the data matrixes do not conform");
   
    nz =
    nx = INTEGER(GET_DIM(x))[0];
    ny = INTEGER(GET_DIM(y))[0];

    if (m == 2 && nx != ny)
	error("the number of rows of the data matrixes do not conform");

    if (TYPEOF(x) != REALSXP) {
	PROTECT(x = coerceVector(R_x, REALSXP));

	if (isNull(R_y) || R_x == R_y)
	    y = x;
    }

    if (TYPEOF(y) != REALSXP) 
	PROTECT(y = coerceVector(R_y, REALSXP));

    if (m == 0) {
	SEXP d;

	l = nx;	
        PROTECT(r = allocVector(REALSXP, l*(l-1)/2));
	setAttrib(r, install("Size"), PROTECT(ScalarInteger(nx)));
	UNPROTECT(1);
	
	if (!isNull(d = getAttrib(x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, PROTECT(mkString("dist")));
	UNPROTECT(1);
    } else
    if (m == 1) {
	SEXP d1, d2;
	
        PROTECT(r = allocMatrix(REALSXP, nx, ny));

	d1 = getAttrib(x, R_DimNamesSymbol);
	d2 = getAttrib(y, R_DimNamesSymbol);
	if (!isNull(d1) || !isNull(d2)) {
	    SEXP d;

	    setAttrib(r, R_DimNamesSymbol, PROTECT(d = allocVector(VECSXP, 2)));
	    UNPROTECT(1);
	    SET_VECTOR_ELT(d, 0, isNull(d1) ? d1 : VECTOR_ELT(d1, 0));
	    SET_VECTOR_ELT(d, 1, isNull(d2) ? d2 : VECTOR_ELT(d2, 0));
	}
    } else
	PROTECT(r = allocVector(REALSXP, nx));
    
    PROTECT(tx = allocVector(REALSXP, n)); 
    PROTECT(ty = allocVector(REALSXP, n));

    PROTECT(c = LCONS(f, CONS(tx, CONS(ty, p))));
    
    l = 0;
    for (j = 0; j < ny; j++) {
        for (k = 0; k < n; k++)
            REAL(ty)[k] = REAL(y)[j+k*ny];
	if (m == 0)
	    i = j+1;
	else 
	if (m == 1)
	    i = 0;
	else {
	    i  = j;
	    nz = j+1;
	}
        for (; i < nz; i++) {
            for (k = 0; k < n; k++)
                REAL(tx)[k] = REAL(x)[i+k*nx];
            SEXP s = eval(c, R_GlobalEnv);
            if (LENGTH(s) != 1)
                error("not a scalar return value");
	    // fixme: warning?
	    if (TYPEOF(s) != REALSXP) {
                REAL(r)[l++] = REAL(coerceVector(PROTECT(s), REALSXP))[0];
		UNPROTECT(1);
	    }
	    else
		REAL(r)[l++] = REAL(s)[0];
        }                       
        R_CheckUserInterrupt();
    }

    UNPROTECT(4);
    if (x != R_x)
	UNPROTECT(1);
    if (!isNull(R_y) && R_y != R_x && y != R_y)
	UNPROTECT(1);

    return r;
}

/* compute auto- or cross-distances with a user-supplied
 * function given list data.
 * 
 * ceeboo 2006, 2007
 */

SEXP R_apply_dist_list(SEXP p) {
    int i, j, nx, ny, nz, m = 0;
    SEXP r, c, d, tx, ty;
    SEXP x, y, f;
    R_xlen_t l;

    p = CDR(p);
    if (length(p) < 4)
        error("invalid number of arguments");
    x = CAR(p); y = CADR(p);
    if (TYPEOF(x) != VECSXP || (!isNull(y) && TYPEOF(y) != VECSXP))
        error("invalid data parameter(s)");
    p = CDDR(p); d = CAR(p);
    if (TYPEOF(d) != LGLSXP)
	error("invalid option parameter");
    p = CDR(p); f = CAR(p); 
    if (!isFunction(f))
        error("invalid function parameter");
    p = CDR(p);

    if (isNull(y))
        y = x;  
    else
    if (LOGICAL(d)[0] == TRUE)
	m = 2;
    else
	m = 1;
  
    nz =
    nx = LENGTH(x);
    ny = LENGTH(y);

    if (m == 0) {
	SEXP d;

	l = nx;	
        PROTECT(r = allocVector(REALSXP, l*(l-1)/2));
	
	setAttrib(r, install("Size"), PROTECT(ScalarInteger(nx)));
	UNPROTECT(1);
	
	PROTECT(d = getAttrib(x, R_NamesSymbol));
	if (!isNull(d))
	    setAttrib(r, install("Labels"), d);
	UNPROTECT(1);
	// fixme: package?
	setAttrib(r, R_ClassSymbol, PROTECT(mkString("dist")));
	UNPROTECT(1);
    } else 
    if (m == 1) {
	SEXP d1, d2;
	
        PROTECT(r = allocMatrix(REALSXP, nx, ny));
	
	PROTECT(d1 = getAttrib(x, R_NamesSymbol));
	PROTECT(d2 = getAttrib(y, R_NamesSymbol));
	if (!isNull(d1) || !isNull(d2)) {
	    SEXP d;

	    setAttrib(r, R_DimNamesSymbol, PROTECT(d = allocVector(VECSXP, 2)));
	    UNPROTECT(1);
	    SET_VECTOR_ELT(d, 0, d1);
	    SET_VECTOR_ELT(d, 1, d2);
	}
	UNPROTECT(2);
    } else {
	if (nx != ny)
	    error("the number of components of 'x' and 'y' does not conform");

	PROTECT(r = allocVector(REALSXP, nx));
    }

    PROTECT(ty = CONS(R_NilValue, p));
    PROTECT(tx = CONS(R_NilValue, ty));
    PROTECT(c = LCONS(f, tx));
   
    l = 0;
    for (j = 0; j < ny; j++) {
        SETCAR(ty, VECTOR_ELT(y, j));
	if (m == 0)
	    i = j+1;
	else
	if (m == 1)
	    i = 0;
	else {
	    i  = j;
	    nz = j+1;
	}
        for (; i < nz; i++) {
            SETCAR(tx, VECTOR_ELT(x, i));
            SEXP s = eval(c, R_GlobalEnv);
            if (LENGTH(s) != 1)
                error("not a scalar return value");
	    // fixme: warning?
	    if (TYPEOF(s) != REALSXP) {
                REAL(r)[l++] = REAL(coerceVector(PROTECT(s), REALSXP))[0];
		UNPROTECT(1);
	    }
	    else
		REAL(r)[l++] = REAL(s)[0];
        }                       
        R_CheckUserInterrupt();
    }

    UNPROTECT(4);

    return r;
}

/* compute binary auto- or cross-distances with a user-supplied
 * function given logical matrix data, and by precomputing the
 * number of concordant and discordant pairs.
 *
 * dm 2007
 */

SEXP R_apply_dist_binary_matrix(SEXP p) {
    int i, j, k, n, nx, ny, nz, m = 0;   
    int i0, j0;
    SEXP r, c, d, ta, tb, tc, td, tn;
    SEXP x, y, f;
    R_xlen_t l;

    p = CDR(p);
    if (length(p) < 3)
        error("invalid number of arguments");
    x = CAR(p); y = CADR(p);
    if (!isMatrix(x) || TYPEOF(x) != LGLSXP || 
	(!isNull(y) && (!isMatrix(y) || TYPEOF(x) != LGLSXP)))
        error("invalid data parameter(s)");
    p = CDDR(p); d = CAR(p);
    if (TYPEOF(d) != LGLSXP)
	error("invalid option parameter");
    p = CDR(p); f = CAR(p); 
    if (!isFunction(f))
        error("invalid function parameter");
    p = CDR(p);

    if (isNull(y))
        y = x;  
    else
    if (LOGICAL(d)[0] == TRUE)
	m = 2;
    else
        m = 1;
    
    if ((n = INTEGER(GET_DIM(x))[1]) != INTEGER(GET_DIM(y))[1])
        error("data parameters do not conform");

    nz =
    nx = INTEGER(GET_DIM(x))[0];
    ny = INTEGER(GET_DIM(y))[0];

    if (m == 0) {
	SEXP d;

	l = nx;	
        PROTECT(r = allocVector(REALSXP, l*(l-1)/2));
	
	setAttrib(r, install("Size"), PROTECT(ScalarInteger(nx)));
	UNPROTECT(1);
	
	if (!isNull(d = getAttrib(x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, PROTECT(mkString("dist")));
	UNPROTECT(1);
    } else
    if (m == 1) {
	SEXP d1, d2;
	
        PROTECT(r = allocMatrix(REALSXP, nx, ny));
	
	d1 = getAttrib(x, R_DimNamesSymbol);
	d2 = getAttrib(y, R_DimNamesSymbol);
	if (!isNull(d1) || !isNull(d2)) {
	    SEXP d;

	    setAttrib(r, R_DimNamesSymbol, PROTECT(d = allocVector(VECSXP, 2)));
	    UNPROTECT(1);
	    SET_VECTOR_ELT(d, 0, isNull(d1) ? d1 : VECTOR_ELT(d1, 0));
	    SET_VECTOR_ELT(d, 1, isNull(d2) ? d2 : VECTOR_ELT(d2, 0));
	}
    } else {
	if (nx != ny)
	    error("the number of rows of 'x' and 'y' does not conform");

	PROTECT(r = allocVector(REALSXP, nx));
    }

    PROTECT(ta = allocVector(INTSXP, 1)); 
    PROTECT(tb = allocVector(INTSXP, 1)); 
    PROTECT(tc = allocVector(INTSXP, 1)); 
    PROTECT(td = allocVector(INTSXP, 1)); 
    PROTECT(tn = allocVector(INTSXP, 1)); 

    PROTECT(c = LCONS(f, CONS(ta, CONS(tb, CONS(tc, CONS(td, CONS(tn, p))))
)));
    
    l = 0;
    for (j = 0; j < ny; j++) {
	if (m == 0)
	    i = j+1;
	else
	if (m == 1)
	    i = 0;
	else {
	    i  = j;
	    nz = j+1;
	}
        for (; i < nz; i++) {
            INTEGER(ta)[0] = INTEGER(tb)[0] = INTEGER(tc)[0] = INTEGER(tn)[0] = 
0;
            for (k = 0; k < n; k++) {
                i0 = LOGICAL(x)[i + k * nx];
                j0 = LOGICAL(y)[j + k * ny];
                if (i0 == NA_LOGICAL || j0 == NA_LOGICAL)
                    continue;
                INTEGER(ta)[0] += (i0 == TRUE  && j0 == TRUE);
                INTEGER(tb)[0] += (i0 == TRUE  && j0 == FALSE);
                INTEGER(tc)[0] += (i0 == FALSE && j0 == TRUE);
                INTEGER(tn)[0]++;
            }
            if (INTEGER(tn)[0] == 0)
                INTEGER(td)[0] = 0;
            else
                INTEGER(td)[0] = INTEGER(tn)[0] - INTEGER(ta)[0] - 
                                                  INTEGER(tb)[0] - 
                                                  INTEGER(tc)[0];
            SEXP s = eval(c, R_GlobalEnv);
            if (LENGTH(s) != 1)
                error("not a scalar return value");
            // fixme: warning?
            if (TYPEOF(s) != REALSXP) {
                REAL(r)[l++] = REAL(coerceVector(PROTECT(s), REALSXP))[0];
		UNPROTECT(1);
	    }
            else
                REAL(r)[l++] = REAL(s)[0];
        }                       
        R_CheckUserInterrupt();
    }
    
    UNPROTECT(7);

    return r;
}

/* compute auto- or cross-distances with a user-supplied
 * function given data.frame data.
 *
 * because of the details this is insane ...
 *
 * ceeboo 2007
 */

static void setElement(SEXP x, int i, SEXP y) {
    switch (TYPEOF(x)) {
    case LGLSXP:
	LOGICAL(x)[0] = LOGICAL(y)[i];
	break;
    case INTSXP:
	INTEGER(x)[0] = INTEGER(y)[i];
	break;
    case REALSXP:
	REAL(x)[0]    = REAL(y)[i];
	break;
    case STRSXP:
	SET_STRING_ELT(x, 0, STRING_ELT(y, i));
	break;
    case VECSXP:
	SET_VECTOR_ELT(x, 0, VECTOR_ELT(y, i));
	break;
    default:
	error("type not implemented");
    }
}

SEXP R_apply_dist_data_frame(SEXP p) {
    int i, j, k, nc, nx, ny, nz, m = 0;
    SEXP r, c, d, tx, ty, rx, ry;
    SEXP x, y, f;
    R_xlen_t l;

    p = CDR(p);
    if (length(p) < 4)
        error("invalid number of arguments");
    x = CAR(p); y = CADR(p);
    if (!inherits(x, "data.frame") || (!isNull(y) && 
        !inherits(y, "data.frame")))
        error("invalid data parameter(s)");
    p = CDDR(p); d = CAR(p);
    if (TYPEOF(d) != LGLSXP)
	error("invalid option parameter");
    p = CDR(p); f = CAR(p); 
    if (!isFunction(f))
        error("invalid function parameter");
    p = CDR(p);

    nc = LENGTH(x);
    if (nc == 0)
	error("cannot handle empty data frames");

    nx = ny = nz = LENGTH(VECTOR_ELT(x, 0));

    if (isNull(y))
        y = x;  
    else {
	if (LENGTH(y) != nc)
	    error("data parameters do not conform");

	ny = LENGTH(VECTOR_ELT(y, 0));

	for (k = 0; k < nc; k++) {
	    if (TYPEOF(VECTOR_ELT(x, k)) != TYPEOF(VECTOR_ELT(y, k)))
		error("data parameters do not conform");
	    // sucks: the c code in identical.c is not
	    //        accessible.
	    c =	eval(PROTECT(LCONS(install("identical"), 
		     PROTECT( CONS(PROTECT(eval(PROTECT(LCONS(install("attributes"), PROTECT(CONS(VECTOR_ELT(x, k), R_NilValue)))), R_GlobalEnv)), 
			      CONS(PROTECT(eval(PROTECT(LCONS(install("attributes"), PROTECT(CONS(VECTOR_ELT(y, k), R_NilValue)))), R_GlobalEnv)),
				   R_NilValue))))), R_GlobalEnv);
	    UNPROTECT(8);
	    if (LOGICAL(c)[0] == FALSE)
		error("attributes of data parameters do not conform");
	}
	if (LOGICAL(d)[0] == TRUE) {
	    if (nx != ny)
		error("the number of rows of 'x' and 'y' do not conform");

	    m = 2;
	} else
	    m = 1;
    }

    // fixme: row.names
    if (m == 0) {
	l = nx;
        PROTECT(r = allocVector(REALSXP, l*(l-1)/2));

	setAttrib(r, install("Size"), PROTECT(ScalarInteger(nx)));
	UNPROTECT(1);
	setAttrib(r, install("Labels"), PROTECT(coerceVector(PROTECT(getAttrib(x, install("row.names"))), STRSXP)));
	UNPROTECT(2);

	setAttrib(r, R_ClassSymbol, PROTECT(mkString("dist")));
	UNPROTECT(1);
    } else 
    if (m == 1) {
	SEXP d;

        PROTECT(r = allocMatrix(REALSXP, nx, ny));

	setAttrib(r, R_DimNamesSymbol, PROTECT(d = allocVector(VECSXP, 2)));
	UNPROTECT(1);
	SET_VECTOR_ELT(d, 0, coerceVector(PROTECT(getAttrib(x, install("row.names"))), STRSXP));
	UNPROTECT(1);
	SET_VECTOR_ELT(d, 1, coerceVector(PROTECT(getAttrib(y, install("row.names"))), STRSXP));
	UNPROTECT(1);
    } else
	PROTECT(r = allocVector(REALSXP, nx));

    PROTECT(tx = allocVector(VECSXP, nc));
    setAttrib(tx, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    setAttrib(tx, install("row.names"), PROTECT(rx = allocVector(INTSXP, 1)));
    UNPROTECT(1);
    setAttrib(tx, R_ClassSymbol, getAttrib(x, R_ClassSymbol));

    PROTECT(ty = allocVector(VECSXP, nc));
    setAttrib(ty, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    setAttrib(ty, install("row.names"), PROTECT(ry = allocVector(INTSXP, 1)));
    UNPROTECT(1);
    setAttrib(ty, R_ClassSymbol, getAttrib(x, R_ClassSymbol));

    for (k = 0; k < nc; k++) {
	SEXP t, s;

	t = VECTOR_ELT(x, k);

					// fixme: should fail for S4
	SET_VECTOR_ELT(tx, k, (s = allocVector(TYPEOF(t), 1)));
	SHALLOW_DUPLICATE_ATTRIB(s, t);	// fixme: may be wrong

	SET_VECTOR_ELT(ty, k, (s = allocVector(TYPEOF(t), 1)));
	SHALLOW_DUPLICATE_ATTRIB(s, t);
    }

    PROTECT(c = LCONS(f, CONS(tx, CONS(ty, p))));

    l = 0;
    for (j = 0; j < ny; j++) {
	for (k = 0; k < nc; k++)
	    setElement(VECTOR_ELT(ty, k), j, VECTOR_ELT(y, k));
	INTEGER(ry)[0] = j+1;		// R index
	if (m == 0)
	    i = j+1;
	else
	if (m == 1)
	    i = 0;
	else {
	    i  = j;
	    nz = j+1;
	}
        for (; i < nz; i++) {
	    for (k = 0; k < nc; k++)
		setElement(VECTOR_ELT(tx, k), i, VECTOR_ELT(x, k));
	    INTEGER(rx)[0] = i+1;
            SEXP s = eval(c, R_GlobalEnv);
            if (LENGTH(s) != 1)
                error("not a scalar return value");
	    // fixme: warning?
	    if (TYPEOF(s) != REALSXP) {
                REAL(r)[l++] = REAL(coerceVector(PROTECT(s), REALSXP))[0];
		UNPROTECT(1);
	    }
	    else
		REAL(r)[l++] = REAL(s)[0];
        }                       
        R_CheckUserInterrupt();
    }

    UNPROTECT(4);

    return r;
}

//
