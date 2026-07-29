#include "logsumexp2.h"


/**
 * @title Exp-Log-Trick for two vectors of the same length
 * @param a n vector
 * @param b n vector
 * @param n length of a,b
 * @return n vector; log(exp(a) + exp(b)) using exp-log-trick
 * @author Erik Hintz
 */
void logsumexp2_c(double *a, double *b, int n, double *res)
{
    int i = 0;
    /* go through input vector */
    for(i = 0; i < n; i++){
        if(a[i] > b[i]){
            /* a[i] is max */
            res[i] = a[i] + log1p(exp(b[i] - a[i]));
        } else {
            /* b[i] is max */
            res[i] = b[i] + log1p(exp(a[i] - b[i]));
        }
    }
}


/**
 * @title R Interface for logsumexp2_c_c()
 * @param see logsumexp2_c() above
 * @return see elogsumexp2_c() above
 * @author Erik Hintz
 */


SEXP logsumexp2(SEXP a, SEXP b, SEXP n)
{
    int n_ = asInteger(n); /* for allocation */
    SEXP res = PROTECT(allocVector(REALSXP, n_)); /* allocate memory*/
    double *res_ = REAL(res); /* pointer to values of res */
    /* Main */
    logsumexp2_c(REAL(a), REAL(b), INTEGER(n)[0], res_);
    UNPROTECT(1);
    return res;
}
