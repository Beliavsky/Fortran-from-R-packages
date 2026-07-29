#include "VaR_bounds.h"


/* The following is working but a) not faster and b) creates problems under
 * Windows if R < 3.3.0 */
#ifdef NOT_NOW_
/**
 * @title Fast order(order(x, decreasing=TRUE))
 * @param x N-vector (one column in the (A)RA()/rearrange() input matrix X)
 * @return order(order(x, decreasing=TRUE))
 * @author Marius Hofert, Kurt Hornik
 */
SEXP indices_opp_ordered_to(SEXP x)
{
    /* Setup */
    int N = LENGTH(x);
    int i;
    SEXP ind = PROTECT(allocVector(INTSXP, N)); /* vector of indices */
    int *ind_; /* define pointer to ind */
    ind_ = INTEGER(ind); /* set pointer to ind */
    SEXP res = PROTECT(allocVector(INTSXP, N)); /* result vector */
    int *res_; /* define pointer to res */
    res_ = INTEGER(res); /* set pointer to res */

    /* Compute order(order(x, decreasing=TRUE)) */
    R_orderVector1(ind_, /* result */
                   N, /* length */
    	           x, /* argument */
    	           TRUE, /* nalast as in order() */
    	           TRUE); /* decreasing TRUE */
    R_orderVector1(res_, /* result */
    	           N, /* length */
    	           ind, /* argument */
    	           TRUE, /* nalast as in order() */
    	           FALSE); /* decreasing FALSE */
    for(i=0; i<N; i++) res_[i] += 1; /* increase all by 1 */

    /* Return */
    UNPROTECT(2);
    return(res);
}
#endif

/**
 * @title Fast split(x, col(x))
 * @param x (N,d)-matrix ((A)RA()/rearrange() input matrix X)
 * @return split(x, col(x))
 * @author Marius Hofert, Kurt Hornik
 */
SEXP col_split(SEXP x)
{
    /* Setup */
    int *dims = INTEGER(getAttrib(x, R_DimSymbol));
    int n = dims[0], d = dims[1];
    SEXP res = PROTECT(allocVector(VECSXP, d));
    SEXP ref;
    int i = 0, j, k; /* i runs globally, j runs over all cols, k runs over all rows */

    /* Distinguish int/real matrices */
    switch (TYPEOF(x)) {
    case INTSXP:
    	for(j = 0; j < d; j++) {
    		SET_VECTOR_ELT(res, j, allocVector(INTSXP, n));
    		int *e = INTEGER(VECTOR_ELT(res, j));
    		for(k = 0 ; k < n ; i++, k++) {
    			e[k] = INTEGER(x)[i];
    		}
    	}
    	break;
    case REALSXP:
    	for(j = 0; j < d; j++) {
    		SET_VECTOR_ELT(res, j, allocVector(REALSXP, n));
    		double *e = REAL(VECTOR_ELT(res, j));
    		for(k = 0 ; k < n ; i++, k++) {
    			e[k] = REAL(x)[i];
    		}
    	}
    	break;
    case LGLSXP:
    	for(j = 0; j < d; j++) {
    		SET_VECTOR_ELT(res, j, allocVector(LGLSXP, n));
    		int *e = LOGICAL(VECTOR_ELT(res, j));
    		for(k = 0 ; k < n ; i++, k++) {
    			e[k] = LOGICAL(x)[i];
    		}
    	}
    	break;
    case STRSXP:
    	for(j = 0; j < d; j++) {
		ref = allocVector(STRSXP, n); /* needs to be here (otherwise the last column is replicated) */
    		SET_VECTOR_ELT(res, j, ref);
    		ref = VECTOR_ELT(res, j);
    		for(k = 0 ; k < n ; i++, k++) {
    			SET_STRING_ELT(ref, k, STRING_ELT(x, i));
    		}
    	}
    	break;
    default: error("Wrong type of 'x': %s", CHAR(type2str_nowarn(TYPEOF(x))));
    }

    /* Return */
    UNPROTECT(1);
    return(res);
}
