#include "precond.h"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

void precond(double *lower, double *upper, double *scale,
             double *C, double *meansqrtmix, double *tol,
             int *d_, int *perm, int *status)

{
    /* Note: element (row, col) for row >= col in {0,...,d-1} has
     index col * (d-1) - col*(col-1)/2 + row */
    double tmp, tmpa, tmpb, sum, denom, scprod, currexplength; /* to store a temporary values (double) */
    int tmpint; /* to store a temporary value, typically index (int) */
    int minind = 0; 
    int d = *d_; /* dereference for readability */
    double *y;
    y = (double *) malloc(d*sizeof(double)); /* to store conditional expected values */

    /* Go through all variables except last */
    for(int j = 0; j < (d-1); j++){

        if(j == 0){
            /* Find argmin_l { <expected length of interval l> } */
            currexplength = 2; /* length always in (0, 1) => initialize to >1 */
            for(int l = 0; l < d; l++){
                /* Compute expected length for variable l */
                denom = sqrt( scale[l*d - l*(l-1)/2]) * meansqrtmix[l];
                tmp = pnorm( upper[l] / denom, 0, 1, 1, 0) - pnorm( lower[l] / denom, 0, 1, 1, 0);

                if(tmp <= currexplength){ /* found index with smaller expected length */
                    currexplength = tmp;
                    minind = l;
                }
            }
            /* Swap variable 0 and 'minind' */
            if(minind > 0){
                swap(lower, upper, meansqrtmix, perm, scale, d, 0, minind);
            }
            /* Update 'C' */
            C[0] = sqrt(scale[0]);
            for(int l = 1; l < d; l++){
                C[l] = scale[l] / C[0];
            }
            /* Store conditional expected value */
            tmpa = lower[0] / (meansqrtmix[0] * C[0]);
            tmpb = upper[0] / (meansqrtmix[0] * C[0]);
            denom = pnorm(tmpb, 0, 1, 1, 0) - pnorm(tmpa, 0, 1, 1, 0);
            if(denom < *tol){  /* Avoid division by ~0 */
                denom = *tol;
            }
            tmp = (dnorm(tmpa, 0, 1, 0) - dnorm(tmpb, 0, 1, 0)) / denom;
            *(y+0) = (dnorm(tmpa, 0, 1, 0) - dnorm(tmpb, 0, 1, 0)) / denom;
        } else {
            /* Case j > 0 */

            /* Find argmin_l { <expected length of interval l> } */
            currexplength = 2; /* length always in (0, 1) => initialize to >1 */
            for(int l = j; l < d; l++){
                /* Compute the scalarproduct sum_{i=0}^{j-1} C[l,i]*y[i]
                   and sum_{i=0}^{j-1} C[l, i]^2 */
                scprod = 0;
                sum = 0;
                for(int i = 0; i < j; i++){
                    tmpint = i*(d-1)-i*(i-1)/2+l;
                    scprod += *(y+i) * C[tmpint];
                    sum    += C[tmpint] * C[tmpint];
                }
                /* Denominator = sqrt(scale[l,l] - sum_{i=0}^{j-1} C[l, i]^2 ) */
                tmp = scale[l*d - l*(l-1)/2] - sum;
                if(tmp <= 0){ /* Can't take root => return and set status = 10+j (error) */
                    *status = 10+j;
                    return;
                } else {
                    denom = sqrt(tmp);
                }
                tmp = pnorm((upper[l]/meansqrtmix[l] - scprod)/denom, 0, 1, 1, 0) -
                         pnorm((lower[l]/meansqrtmix[l] - scprod)/denom, 0, 1, 1, 0);
                if(tmp < currexplength){ /* found index with smaller expected length */
                    currexplength = tmp;
                    minind = l;
                }
            }

            if(j != minind){
                /* Swap variable 'j' and 'minind'  */
                swap(lower, upper, meansqrtmix, perm, scale, d, j, minind);
                /* Swap row 'j' and 'minind' >= j in the cholesky factor C */
                /* Only the first 'j' columns need to be considered (lower triangle)*/
                for(int l = 0; l < j; l++){
                    /* Swap C[minind, l] and C[j, l]*/
                    tmp = C[l*(d-1)-l*(l-1)/2+minind];
                    C[l*(d-1)-l*(l-1)/2+minind] = C[l*(d-1)-l*(l-1)/2+j];
                    C[l*(d-1)-l*(l-1)/2+j] = tmp;
                }
            }
            /* Update cholesky factor 'C' */
            sum = 0;
            /* Compute sum_{i=0}^{j-1} C[j,i]^2 */
            for(int i = 0; i < j; i++){
                tmpint = i*(d-1)-i*(i-1)/2+j;
                sum += C[tmpint] * C[tmpint];
            }
            /*Rprintf("Value of sum: %f\n", sum);*/
            /* Update C[j,j] */
            tmpint = j*d-j*(j-1)/2;
            C[tmpint] = sqrt(scale[tmpint] - sum);
            /* Update C[l,j] for l > j */
            for(int l = j+1; l < d; l++){
                /* Compute sum_{i=0}^{j-1} C[j,i] C[l,i] */
                sum = 0;
                for(int i = 0; i < j; i++){
                    sum += C[i*(d-1)-i*(i-1)/2+j] * C[i*(d-1)-i*(i-1)/2+l];
                }
                C[j*(d-1)-j*(j-1)/2+l] = (scale[j*(d-1)-j*(j-1)/2+l] - sum)/C[j*d-j*(j-1)/2];
            }
            /* Store conditional expected value */
            /* Compute scalarproduct sum_{l=0}^{j-1} C[j,l] y[l] */
            scprod = 0;
            for(int l = 0; l < j; l++){
                scprod += y[l] * C[l*(d-1)-l*(l-1)/2+j];
            }
            tmpa = (lower[j]/meansqrtmix[j] - scprod)/C[j*d-j*(j-1)/2];
            tmpb = (upper[j]/meansqrtmix[j] - scprod)/C[j*d-j*(j-1)/2];
            denom = pnorm(tmpb, 0, 1, 1, 0) - pnorm(tmpa, 0, 1, 1, 0);
            if(denom < *tol){ /* Avoid division by ~0 */
                denom = *tol;
            }
            *(y+j)= (dnorm(tmpa, 0, 1, 0) - dnorm(tmpb, 0, 1, 0)) / denom;
            
        }
    } /* for (j ..) */
    /* Update C[d,d] */
    /* Compute sum_{i=0}^{d-2} C[d, i]^2 */
    free(y);
    sum = 0;
    for(int i = 0; i < d-1; i++){
        tmpint = (i+1)*(d-1)-i*(i-1)/2;
        sum += C[tmpint] * C[tmpint];
    }
    tmp = scale[d*(d+1)/2 - 1] - sum;
    if(tmp > *tol){
        C[d*(d+1)/2 - 1] = sqrt(tmp);
    } else {
        *status = 2; /* recompute chol(scale) in R */
    }
}



void swap(double *lower, double *upper, double *meansqrtmix, int *perm, double *scale,
           int d, int i, int j)
{
    
    /* To store temporary values for swapping */
    double tmp;
    int tmpint;
    /* Set wlog i < j (swapping i,j <=> swapping j,i) */
    if(j < i){
        tmpint = i;
        i = j;
        j = tmpint;
    }

    /* Swap elements i and j in 'lower', 'upper', 'perm' and 'meansqrtmix' */
    tmp = lower[i];
    lower[i] = lower[j];
    lower[j] = tmp;
    
    tmp = upper[i];
    upper[i] = upper[j];
    upper[j] = tmp;
    
    tmpint = perm[i];
    perm[i] = perm[j];
    perm[j] = tmpint;
    
    tmp = meansqrtmix[i];
    meansqrtmix[i] = meansqrtmix[j];
    meansqrtmix[j] = tmp;
    

    /* Swap cov(Xi, Xk) <-> Cov(Xj, Xk) for all k */
    for(int k = 0; k < d; k++){
        if(k < i){
            /* Swap (i,k) and (j,k) */
            tmp = scale[k*(d-1)-k*(k-1)/2+i];
            scale[k*(d-1)-k*(k-1)/2+i] = scale[k*(d-1)-k*(k-1)/2+j];
            scale[k*(d-1)-k*(k-1)/2+j] = tmp;
        } else if(k == i){ /* Var(Xi) <-> Var(Xj) once! */
            /* Swap (i,i) and (j,j) */
            tmp = scale[i*d-i*(i-1)/2];
            scale[i*d-i*(i-1)/2] = scale[j*d-j*(j-1)/2];
            scale[j*d-j*(j-1)/2] = tmp;
        } else if(i < k && k < j){
            /* Swap (k,i) and (j,k) */
            tmp = scale[i*(d-1)-i*(i-1)/2+k];
            scale[i*(d-1)-i*(i-1)/2+k] = scale[k*(d-1)-k*(k-1)/2+j];
            scale[k*(d-1)-k*(k-1)/2+j] = tmp;
        } else if(k > j){
            /* Swap (k,i) and (k,j) */
            tmp = scale[i*(d-1)-i*(i-1)/2+k];
            scale[i*(d-1)-i*(i-1)/2+k] = scale[j*(d-1)-j*(j-1)/2+k];
            scale[j*(d-1)-j*(j-1)/2+k] = tmp;
        }
    }
    
}
