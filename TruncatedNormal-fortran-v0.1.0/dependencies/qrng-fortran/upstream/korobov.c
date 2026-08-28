/* C function for computing a Korobov sequence ********************************/

#include "korobov.h"


/**
 * @title Generate n Points of a d-dimensional Korobov Sequence
 * @param n number of points
 * @param d dimension
 * @param generator vector of generator points
 * @param randomize string indicating whether the points are randomized
 * @param res pointer to the result matrix
 * @return void
 * @author Marius Hofert based on C. Lemieux's RandQMC
 */
void korobov(int n, int d, int *generator, int randomize, double *res)
{
	int i, j, ij;
	double U;
	double *aux;
	aux = (double *) R_alloc(d, sizeof(double));

	/* Init */
	for(j=0; j<d; j++){
		aux[j] = generator[j] / ((double) n);
		res[j*n] = 0.0; /* case i = 0 below */
	}

	/* Generate points */
	for(i=1; i<n; i++){ /* omit i=0 as done in init above */
		for(j=0; j<d; j++){
			ij = j*n+i;
			res[ij] = res[j*n + (i-1)] + aux[j];
			if(res[ij] > 1) res[ij] = res[ij] - 1.0;
		}
	}

	/* Randomization */
	if(randomize == 1) {
		GetRNGstate();
		for(j=0; j<d; j++){
			U = unif_rand();
			for(i=0; i<n; i++){
				ij = j*n+i;
				res[ij] = res[ij] + U;
				if(res[ij] > 1) res[ij] = res[ij] - 1.0;
			}
		}
		PutRNGstate();
	}
}

/**
 * @title R Interface to C for Generating a Korobov Sequence
 * @param n_ number of points
 * @param d_ dimension
 * @param generator_ vector of generator points
 * @param randomize_ string indicating whether the points are randomized
 * @return (n, d)-matrix
 * @author Marius Hofert
 */
SEXP korobov_(SEXP n_, SEXP d_, SEXP generator_, SEXP randomize_)
{
    /* Input parameters */
    int n = asInteger(n_); /* numeric(1) */
    int d = asInteger(d_); /* numeric(1) */
    int *generator = INTEGER(coerceVector(generator_, INTSXP)); /* numeric(d) */
    int randomize = asLogical(randomize_); /* numeric(1) */

    /* Create result object */
    SEXP res_ = PROTECT(allocMatrix(REALSXP, n, d)); /* (n,d)-matrix */
    double *res = REAL(res_); /* pointer to the values of res */

    /* Main */
    korobov(n, d, generator, randomize, res);

    /* Return */
    UNPROTECT(1); /* clean-up */
    return res_;
}
