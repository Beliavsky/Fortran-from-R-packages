#include "eval_nvmix_integral.h"


/**
 * @title Evaluate Integrand for a Normal Variance Mixture
 * @param lower d-vector of lower evaluation limits
 * @param upper d-vector of upper evaluation limits
 * @param groupings d-vector giving groupings (i.e., component i has mixing rv groupings[i])
 * @param numgroups integer, number of groups (=1 for NVM, >1 for gNVM)
 * @param U vector representing '(n, r-1)-matrix' of uniforms (e.g. Sobol pointset)
 * @param rtW  vector representing '(n, numgroups)' matrix of realizations of sqrt(mix)
 * @param rtWant vector representing '(n, numgroups)' matrix of antitthetic realizations of sqrt(mix)
 * @param n sample size (i.e., number of rows of U)
 * @param d dimension of the original problem
 * @param r rank of factor
 * @param kfactor vector of length r giving height of each step in 'factor'
 * @param factor lower triangular Cholesky factor as vector (factor[row, col] = factor[col*(d-1)-col*(col-1)/2+row])
 * @param ZERO smallest number x > 0 such that x != 0
 * @param ONE   largest number x < 1 such that x != 1
 * @param res 2-vector to store results (mean + estimated variance)
 * @return 2-vector consisting of mean(f(U)) and var(f(U)) using antithetic variates
 * @note lower and upper can have -/+Inf entries. While pnorm() would give the
 *       correct result, it safes time to check each time if the argument is -/+Inf
 *       and setting the value to 0/1 rather than calling pnorm().
 *       Access factor[row, col] as factor[col*(d-1)-col*(col-1)/2+row]
 * @author Erik Hintz and Marius Hofert
 */


void eval_nvmix_integral_c(double *lower, double *upper, int *groupings, int numgroups,
                           double *U, double *rtW, double *rtWant,
                           int n, int d, int r, int *kfactor, double *factor,
                           double ZERO, double ONE, int doant, double *res)
{
    double yorg[r-1], dorg, difforg, scprodorg, ldorg, ldifforg, lforg;
    /* Antithetic equivalents: */
    double yant[r-1], ldant, ldiffant;
    double lfant = 0;
    double scprodant = 0;
    double diffant = 0;
    double dant = 0;
    /* Note: <name>org stands for "original", <name>ant for antithetic */
    /* l<name>org/ant stands for log of that variable */
    /* y:       vector to save phi^{-1}(dj+uj(ej-dj)) */
    /* d:       current values of di from the paper */
    /* diff:    current values of (ei-di) from the paper */
    /* f:       current value of (e1-d1) * (e2-d2) * ... * (ei-di) */
    /* scprod:  scalar product sum factor_{ij} y_j */
    
    int tmpint; /* to store temporary values (integer) */
    double tmp; /* to store temporary values (double) */
    double lowermaxorg, upperminorg, scprodorgnew;  /* needed in singular case */
    double lowermaxant = 0;
    double upperminant = 0;
    double scprodantnew = 0;
    double sum = 0; /* to store sum_{i=1}^n (y_i + yant_i)/2 */
    double sumsq = 0; /* to store sum_{i=1}^n (y_i + yant_i)^2/4 (for variance calculation) */
    int current_limit; /* index of current element in 'lower'/'upper' */
    int ind_W; /* current index for W */
    int i, j, l, m; /* counters for loops */
    /* Avoid recalculation */
    double qnormzero = qnorm(ZERO, 0, 1, 1, 0);
    double qnormone  = -qnormzero;
    
    /* For each row of U (for access, use U[s,k] = U[k * numrows + s]) */
    for(j = 0; j < n; j++){
        /* Initialize current_limit */
        current_limit = 0;
        /* Grab 'lower' and 'upper' */
        lowermaxorg = lower[current_limit];
        upperminorg = upper[current_limit];
        /* Non-singular case: */
        if(r == d){
            lowermaxorg = lowermaxorg / factor[0];
            upperminorg = upperminorg / factor[0];
        } else if(kfactor[0] > 1){
            /* Singular case: Find active limit */
            for(i = 1; i <= kfactor[0]; i++){
                if(lower[i] > lowermaxorg)
                    lowermaxorg = lower[i];
                if(upper[i] < upperminorg)
                    upperminorg = upper[i];
            }
        }
        
        ind_W = groupings[current_limit] - 1; /* column index for W */
        
        /* Need essentially log( pnorm(upper*) - pnorm(lower*)) */
        /* For higher accuracy, make use of 'lower.tail' argument in pnorm() when helpful */
        /* Note the arguments: .Call(C_pnorm, q, mean, sd, lower.tail, log.p) */
        
        /* lower_ = -Inf */
        if(lowermaxorg == R_NegInf){
            dorg = 0;
            if(doant){
                dant = 0;
            }
            if(upperminorg == R_PosInf){
                /* Case Phi(Inf) - Phi(-Inf) */
                difforg = 1;
                ldifforg = 0;
                if(doant){
                    diffant = 1;
                    ldiffant = 0;
                }
            } else {
                /* Case Phi(upper) - Phi(-Inf) = Phi(upper) */
                ldifforg = pnorm(upperminorg / rtW[ind_W * n + j], 0, 1, 1, 1);
                difforg = exp(ldifforg);
                if(doant){
                    ldiffant = pnorm(upperminorg / rtWant[ind_W * n + j], 0, 1, 1, 1);
                    diffant = exp(ldiffant);
                }
            }
        } else {
            /* lower_ != - Inf */
            if(upperminorg == R_PosInf){
                /* Case Phi(Inf) - Phi(lower_) = Phi(lower_, lower.tail = FALSE) */
                ldifforg = pnorm(lowermaxorg / rtW[ind_W * n + j], 0, 1, 0, 1);
                difforg  = exp(ldifforg);
                dorg = 1 - difforg;
                if(doant){
                    ldiffant = pnorm(lowermaxorg / rtWant[ind_W * n + j], 0, 1, 0, 1);
                    diffant  = exp(ldiffant);
                    dant = 1 - diffant;
                }
            } else {
                /* Case Phi(upper_) - Phi(lower_)  */
                ldorg = pnorm(lowermaxorg / rtW[ind_W * n + j], 0, 1, 1, 1);
                dorg  = exp(ldorg);
                /* logsumexp trick for log(Phi(upper_)-Phi(lower_))*/
                tmp = pnorm(upperminorg / rtW[ind_W * n + j], 0, 1, 1, 1);
                ldifforg = tmp + log1p(-exp(ldorg - tmp));
                difforg = exp(ldifforg);
                if(doant){
                    ldant = pnorm(lowermaxorg / rtWant[ind_W * n + j], 0, 1, 1, 1);
                    dant  = exp(ldant);
                    tmp = pnorm(upperminorg / rtWant[ind_W * n + j], 0, 1, 1, 1);
                    ldiffant = tmp + log1p(-exp(ldant - tmp));
                    diffant = exp(ldiffant);
                }
            }
        }
        
        current_limit += kfactor[0];
        /* Go through all r-1 columns of U */
        lforg = ldifforg;
        if(doant){
            lfant = ldiffant;
        }
        for(i = 0; i < r-1; i++){
            ind_W = groupings[current_limit] - 1;
            /* U[i * n + j] corresponds to U[j, i] in the orginal matrix */
            tmp = dorg + U[i * n + j] * difforg;
            /* Check if too close to 0 or 1 */
            if(tmp < ZERO){
                yorg[i] = qnormzero;
            } else if(tmp > ONE) {
                yorg[i] = qnormone;
            } else {
                yorg[i] = qnorm(tmp, 0, 1, 1, 0);
            }
            scprodorg = 0;
            
            if(doant){
                /* The same for the antithetic value */
                tmp = dant + (1-U[i * n + j]) * diffant;
                if(tmp < ZERO){
                    yant[i] = qnormzero;
                } else if(tmp > ONE){
                    yant[i] = qnormone;
                } else {
                    yant[i] = qnorm(tmp, 0, 1, 1, 0);
                }
                scprodant = 0;
            }

            /* Calculate the scalar product sum factor[i,j] y[j] for j = 1:(i-1) */
            for(l = 0; l < (i+1); l++){
                tmpint = l*(d-1)-l*(l-1)/2+current_limit;
                scprodorg += yorg[l] * factor[tmpint];
                if(doant){
                    scprodant += yant[l] * factor[tmpint];
                }
            }
            
            lowermaxorg = (lower[current_limit] / rtW[ind_W * n + j] - scprodorg);
            upperminorg = (upper[current_limit] / rtW[ind_W * n + j] - scprodorg);
            if(doant){
                lowermaxant = (lower[current_limit] / rtWant[ind_W * n + j] - scprodant);
                upperminant = (upper[current_limit] / rtWant[ind_W * n + j] - scprodant);
            }

            /* Non-singular case: */
            if(r == d){
                /* Divide by C[i,i] != 0 */
                tmpint = current_limit*d - current_limit*(current_limit-1)/2;
                lowermaxorg = lowermaxorg / factor[tmpint];
                upperminorg = upperminorg / factor[tmpint];
                if(doant){
                    lowermaxant = lowermaxant / factor[tmpint];
                    upperminant = upperminant / factor[tmpint];
                }
            } else if(kfactor[i+1] > 1){
                /* Singular case: Go through the next rows of factor, adjust limit. Note: In the
                 singular case, 'factor's right-most element is always 1. */
                for(m = 1; m < (kfactor[i+1]+1); m++){
                    /* Calculate scalarproduct factor[i+l,j] y[j] for j = 1:(i-1) (next row of 'factor') */
                    scprodorgnew = 0;
                    if(doant){
                        scprodantnew = 0;
                    }
                    for(l = 0; l < (i+1); l++){
                        tmpint = l*(d-1) - l*(l-1)/2 + current_limit + m;
                        scprodorgnew += yorg[l] * factor[tmpint];
                        if(doant){
                            scprodantnew += yant[l] * factor[tmpint];
                        }
                    }
                    /* Update 'lowermaxorg' if necessary */
                    tmp = (lower[current_limit + m] / rtW[(groupings[current_limit+m]-1) * n + j] - scprodorgnew);
                    if(tmp > lowermaxorg){
                        lowermaxorg = tmp;
                    }
                    /* Update 'upperminorg' if necessary */
                    tmp = (upper[current_limit + m] / rtW[(groupings[current_limit+m]-1) * n + j] - scprodorgnew);
                    if(tmp < upperminorg){
                        upperminorg = tmp;
                    }
                    if(doant){
                        /* Update 'lowermaxant' if necessary */
                        tmp = (lower[current_limit + m] / rtWant[(groupings[current_limit+m]-1) * n + j] - scprodantnew);
                        if(tmp > lowermaxant){
                            lowermaxant = tmp;
                        }
                        /* Update 'upperminant' if necessary */
                        tmp = (upper[current_limit + m] / rtWant[(groupings[current_limit+m]-1) * n + j] - scprodantnew);
                        if(tmp < upperminant){
                            upperminant = tmp;
                        }
                    }
                }
            }
            /* Calculate new d */
            /* Note: lower/upper<...>org Inf/-Inf <=> lower/upper<...>ant Inf/-Inf*/
            /* lower_ = -Inf */
            if(lowermaxorg == R_NegInf){
                dorg = 0;
                if(doant){
                    dant = 0;
                }
                if(upperminorg == R_PosInf){
                    /* Case Phi(Inf) - Phi(-Inf) */
                    difforg = 1;
                    ldifforg = 0;
                    if(doant){
                        diffant = 1;
                        ldiffant = 0;
                    }
                } else {
                    /* Case Phi(upper) - Phi(-Inf) = Phi(upper) */
                    ldifforg = pnorm(upperminorg, 0, 1, 1, 1);
                    difforg = exp(ldifforg);
                    if(doant){
                        ldiffant = pnorm(upperminant, 0, 1, 1, 1);
                        diffant = exp(ldiffant);
                    }
                }
            } else {
                /* lower_ != -Inf */
                if(upperminorg == R_PosInf){
                    /* Case Phi(Inf) - Phi(lower_) = Phi(lower_, lower.tail = FALSE) */
                    ldifforg = pnorm(lowermaxorg, 0, 1, 0, 1);
                    difforg  = exp(ldifforg);
                    dorg = 1 - difforg;
                    if(doant){
                        ldiffant = pnorm(lowermaxant, 0, 1, 0, 1);
                        diffant  = exp(ldiffant);
                        dant = 1 - diffant;
                    }
                } else {
                    /* Case Phi(upper_) - Phi(lower_)  */
                    ldorg = pnorm(lowermaxorg, 0, 1, 1, 1);
                    dorg  = exp(ldorg);
                    /* logsumexp trick for log(Phi(upper_)-Phi(lower_))*/
                    tmp = pnorm(upperminorg, 0, 1, 1, 1);
                    ldifforg = tmp + log1p(-exp(ldorg - tmp));
                    difforg = exp(ldifforg);
                    if(doant){
                        ldant = pnorm(lowermaxant, 0, 1, 1, 1);
                        dant  = exp(ldant);
                        tmp = pnorm(upperminant, 0, 1, 1, 1);
                        ldiffant = tmp + log1p(-exp(ldant - tmp));
                        diffant = exp(ldiffant);
                    }
                }
            }
            /* Update products 'forg' and 'fant' */
            lforg += ldifforg;
            if(doant){
                lfant += ldiffant;
            }
            /* Update i in the singular case (as some rows may need to be skipped) */
            current_limit += kfactor[i+1];
        }
        if(doant){
            tmp = (exp(lforg) + exp(lfant))/2;
        } else {
            tmp = exp(lforg);
        }
        sum += tmp;
        sumsq += tmp*tmp;
    }
    res[0] = sum/n;
    res[1] = (sumsq - n * res[0] * res[0])/(n-1);
}

/**
 * @title R Interface for eval_nvmix_integral_c()
 * @param see eval_nvmix_integral_c()
 * @return mean(f(U)) where f is the integrand and U specifies the point-set
 * @author Erik Hintz, Marius Hofert (polishing)
 */
SEXP eval_nvmix_integral(SEXP lower, SEXP upper, SEXP groupings, SEXP numgroups,
                         SEXP U, SEXP rtW, SEXP rtWant, SEXP n, SEXP d, SEXP r,
                         SEXP kfactor, SEXP factor, SEXP ZERO, SEXP ONE, SEXP doant)
{
    SEXP res = PROTECT(allocVector(REALSXP, 2)); /* allocate memory */
    double *res_ = REAL(res); /* pointer to values of res */

    /* Main */
    eval_nvmix_integral_c(REAL(lower), REAL(upper), INTEGER(groupings), INTEGER(numgroups)[0],
                          REAL(U), REAL(rtW), REAL(rtWant), INTEGER(n)[0],
                          INTEGER(d)[0], INTEGER(r)[0], INTEGER(kfactor),
                          REAL(factor), REAL(ZERO)[0], REAL(ONE)[0], INTEGER(doant)[0], res_);
    UNPROTECT(1);
    /* Return */
    return res;
}
