#include "eval_densmix_integrand.h"


/**
 * @title Evaluate Integrand for a the log-density of a normal variance mixture
 * @param W current_n-vector of realizations of W (*sorted*)
 * @param maha2_2 n-vector of mahalabonis distances (maha^2/2) (*sorted*)
 * @param current_n RQMC sample size to estimate the density (length of W)
 * @param n number of evaluation points (length of maha2_2)
 * @param d dimension of the nvmix distribution
 * @param k see formula below; usually d = k
 * @param lconst n-vector of log-constants of the density
 * @param ldensities n-vector to store ldensity values
 * @param c current_n vector to store calculated values
 * @note This function performs the equvivalent of the following R-Code
 *
 *       b <- lconst - k/2 * log(W) - outer(1/W, maha2_2)
 *       bmax <- apply(b, 2, max) # n-vector of maximal b's
 *       ldensities <- - log(current_n) + bmax + log(colSums(exp(b - rep(bmax, each = current_n))))
 *
 * @author Erik Hintz
 */
void eval_densmix_integrand_c(double *W, double *maha2_2, int current_n, int n,
                              int d, int k, double *lconst, double *ldensities,
                              double *c)
{

    /* (Temporary) variables to store current realizations/data point */
    double current_maha, current_lconst, current_c, next_c, current_W, sum_expc;
    double c_max = 0;
    /* Some more counters and indicators */
    int maxindex = 0;
    int found_max;
    int i, j, l;
    int startindex = 0; /* index to start looking for the maximum */
    
    /* Some constants that we can re-use: */
    double neglogcurrent_n = -log( (double) current_n);
    double k2 =  (double) k / 2;
    
    /* For each evaluation point in maha2_2 */
    for(j = 0; j < n; j++){
        /* Grab current maha2_2 */
        current_maha = maha2_2[j];
        current_lconst = lconst[j];
        /* Calculate c_i starting from i = startindex until max reached */
        i = startindex;
        /* Case startindex == current_n -1 special: */
        if(i == current_n - 1){
            /* in this case, maximum found (last row) */
            found_max = 1;
            maxindex = i;
            c_max = current_lconst - k2 * log(W[i]) - current_maha/W[i];
        } else {
            /* if not, need to find the maximum */
            found_max = 0;
            /* first c_i: */
            current_W = W[i];
            current_c = current_lconst - k2 * log(current_W) - current_maha/current_W;
            /* save this so that we can re-use it later. To save memory, we always
             fill the vector c from the beginning, no matter what. */
            c[i - startindex] = current_c;
            i += 1;
            while( (found_max == 0) && (i < current_n)){
                /* Get next realization and safe it */
                current_W = W[i];
                next_c = current_lconst - k2 * log(current_W) - current_maha/current_W;
                c[i - startindex] = next_c;
                
                /* Did we find the maximum? c is first increasing, then decreasing
                 as a function of W. Use "<" so that the case of several W being
                 the same is accounted for (eg two-point-mixture) */
                if(next_c < current_c){
                    found_max = 1;
                    c_max = current_c;
                    /* The maximum occured in the previous index */
                    maxindex = i - 1;
                } else if(i == current_n - 1){
                    /* In case no max was found until here it is the last index: */
                    maxindex = i;
                    found_max = 1;
                    c_max = next_c;
                } else {
                    current_c = next_c;
                    i += 1;
                }
            }
        }
        /* Calculate sum_1^current_n exp(ci-cmax): We omit the i=max term */
        sum_expc = 0;
        for(l = 0; l < current_n; l++){
            /* Was c_l already calculated? */
             if( (l >= startindex) && (l < maxindex ) ){
                sum_expc += exp(c[l - startindex] - c_max);
            } else if( l != maxindex){
                current_W = W[l];
                sum_expc += exp( current_lconst - k2 * log(current_W) - current_maha/current_W - c_max);
            }
            
        }
        /* Position of c_max is increasing in maha2_2. Since the latter is
         sorted in increasing order, we use startindex = maxindex  for the
         next maha2_2 value */
        startindex = maxindex;
        /* Done for maha2_2[j] */
        /* Use log1p(sum_expc) = log(1 + sum_expc) as the exp(0) term was omitted
         when determining sum_expc */
        ldensities[j] = neglogcurrent_n + c_max + log1p(sum_expc);
    }
}


/**
 * @title R Interface for eval_densmix_integrand_c()
 * @param see eval_densmix_integrand_c() above
 * @return see eval_densmix_integrand_c() above
 * @author Erik Hintz
 */


SEXP eval_densmix_integrand(SEXP W, SEXP maha2_2, SEXP current_n, SEXP n, SEXP d, SEXP k, SEXP lconst)
{
    int n_ = asInteger(n); /* for allocation */
    int current_n_ = asInteger(current_n); /* for allocation */
    
    SEXP ldensities = PROTECT(allocVector(REALSXP, n_)); /* allocate memory*/
    double *ldensities_ = REAL(ldensities); /* pointer to values of ldensities */
    
    SEXP c = PROTECT(allocVector(REALSXP, current_n_));
    double *c_ = REAL(c);
    /* Main */
    
    eval_densmix_integrand_c(REAL(W), REAL(maha2_2),
				       INTEGER(current_n)[0], INTEGER(n)[0], INTEGER(d)[0],
				       INTEGER(k)[0], REAL(lconst), ldensities_, c_);
    
    UNPROTECT(2);
    /* Return */
    return ldensities;
}
