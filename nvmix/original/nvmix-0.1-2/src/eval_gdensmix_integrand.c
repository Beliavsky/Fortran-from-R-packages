#include "eval_gdensmix_integrand.h"

/**
 * @title Compute Mahalanobis distances for grouped mixtures
 * @param x (N, d) matrix of density evaluation points
 * @param mix (n, numgroups) matrix of realizations of the mixing variables
 * @param groupings d-vector specifying group-structure
 * @param factorinv lower triangular part of inverse of 'factor' (as vector)
 * @param d integer, dimension
 * @param N integer, number of evaluation points
 * @param n integer, number of mixing realizations
 * @param lconst logarithmic constant to be added (typically -d/2*log(2*pi) - lrdet)
 * @param returnall integer if all logh(u) returned (1) or if LSE-trick is to be applied (0)
 * @return if returnall = 1: (n, N) matrix with log h(u) values; otherwise N-vector with estimated log-density values
 * @author Erik Hintz
 */

/*  Access factor[row, col] as factor[col*(d-1)-col*(col-1)/2+row] */
/* For each row of U (for access, use U[row,col] = U[col * numrows + row]) */
    
void eval_gdenxmis_integrand_c(double *x, double *mix, int *groupings,
            double *factorinv, int d, int N, int n,  double lconst,
            int returnall, double *res)
{
    /* Declarations */
    int i, j, k, l;
    double tmpmaha, tmpsum, lhval;
    double currmax = 0;
    double tmpsumlogmix = 0;
    int currmaxind = 0;
    /* Vectors */
    double *tempx;
    double *templhval;
    double *sumlogmixhalf;
    /* Allocate memory */
    tempx = (double*) malloc(d * sizeof(double)); /* to store a temp value of 'x' */
    templhval = (double*) malloc(n * sizeof(double)); /* to store n values logh(u) for LSE */
    sumlogmixhalf = (double*) malloc(n * sizeof(double)); /* to store sumlogmix for each row in 'mix' */
    

    /* For each row in 'x' */
    for(i = 0; i < N; i++){
        /* For each row in 'mix' */
        for(j = 0; j < n; j++){
            /* Compute maha-distance and sumlogmix */
            tmpmaha = 0;
            if(i == 0){
                tmpsumlogmix = 0;
            }
            for(k = 0; k < d; k++){
                *(tempx+k) = x[k * N + i] / sqrt( mix[(groupings[k]-1) * n + j] );
                /* k'th element in the matrix product factorinv %*% tempx */
                tmpsum = 0;
                for(l = 0; l <= k; l++){
                    tmpsum += factorinv[l*(d-1)-l*(l-1)/2+k] * *(tempx+l);
                }
                tmpmaha += tmpsum*tmpsum;
                if(i == 0){
                    tmpsumlogmix += log(mix[(groupings[k]-1) * n + j]);
                }
            }
            if(i == 0){
                *(sumlogmixhalf+j) = tmpsumlogmix/2;
            }
            lhval = lconst - *(sumlogmixhalf+j) - tmpmaha/2;
            if(returnall == 1){
                res[i*n + j] = lhval;
            } else {
                if(j == 0){
                    currmax = lhval;
                    currmaxind = 0;
                } else {
                    if(lhval > currmax){
                        currmax = lhval;
                        currmaxind = j;
                    }
                }
                *(templhval+j) = lhval;
            } /* else () */
        } /* for(j ..) */
        if(returnall == 0){
            /* Apply log-sum-exp trick */
            /* Compute currmax + log(sum(exp(lhvals[j] - currmax))) */
            tmpsum = 0;
            for(k = 0; k < n; k++){
                if(k != currmaxind){
                    tmpsum += exp(templhval[k] - currmax);
                }
            }
            res[i] = currmax + log1p(tmpsum);
        }
    } /* for(i ..) */
    
    /* Free allocated memory */
    free(tempx);
    free(templhval);
    free(sumlogmixhalf);
    
}


/**
 * @title R Interface for eval_gdenxmis_integrand_c() with 'returnall = 1'
 * @param see above
 * @return (n, N) matrix with values logh(u)
 * @author Erik Hintz
 */
SEXP eval_gdensmix_integrand_returnall(SEXP x, SEXP mix, SEXP groupings, SEXP factorinv,
                                    SEXP d, SEXP N, SEXP n, SEXP lconst)
{
      

    int n_ = asInteger(n);
    int N_ = asInteger(N);
    int returnall = 1;
    /* Allocate memory */
    SEXP res = PROTECT(allocMatrix(REALSXP, n_, N_));
    double *res_ = REAL(res); /* pointer to values of res */

    /* Main */
    eval_gdenxmis_integrand_c(REAL(x), REAL(mix), INTEGER(groupings), REAL(factorinv),
                              INTEGER(d)[0], N_, n_, REAL(lconst)[0], returnall, res_);
    UNPROTECT(1);
    /* Return */
    return res;
}

/**
* @title R Interface for eval_gdenxmis_integrand_c() with 'returnall = 0'
* @param see above
* @return N-vector with estimated log-density values
* @author Erik Hintz
*/
SEXP eval_gdensmix_integrand_LSE(SEXP x, SEXP mix, SEXP groupings, SEXP factorinv,
                                    SEXP d, SEXP N, SEXP n, SEXP lconst)
{
    int n_ = asInteger(n);
    int N_ = asInteger(N);
    int returnall = 0;
    /* Allocate memory */
    SEXP res = PROTECT(allocVector(REALSXP, N_));
    double *res_ = REAL(res); /* pointer to values of res */

    /* Main */
    eval_gdenxmis_integrand_c(REAL(x), REAL(mix), INTEGER(groupings), REAL(factorinv),
                              INTEGER(d)[0], N_, n_, REAL(lconst)[0], returnall,
                              res_);
    UNPROTECT(1);
    /* Return */
    return res;
}



