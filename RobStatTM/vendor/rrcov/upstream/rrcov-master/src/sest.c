#define STRICT_R_HEADERS 1

/* TO DO -  FIXME
 * - Fix the parameters maxit, rtol and converged
 * - Fix "if(ires > 0) ... else"  -there is no need of such distinction
 * - Use QR to compute the inverse of the cov matrix - see mtxdet().
 */
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Applic.h>	/* for the QR	  routines */ 

#define MAX_NTRY 1000 /* was 500 */
#define ZERO 1e-10
#define INFI 1e+20

void r_sample(int *x, int *n, int *k);
void selectwr(int *array, int size, int nchoices);
void reverse(int *a, int n);
void resample(int *array, int n, int k);

void covp(double **x, int *n, int *p, int *indices, int *nind, 
   	  double *mean, double **cov, double *det, int *rank); 
void covpold(double **a, int n, int p, int *id, int np, double *t, double ** cov);
void covar(double **a, int n, int p, double *t, double ** cov);
void covwt(double **a, int n, int p, double *wts, double *t, double **cov);
double mymed(int n, double *x);
double mymedabs(int n, double *x);
double mad(int n, double *x, double *dwork1, double *dwork2, double *mu);
int maxind(double *a, double *maxval, int n);
double norm (double *x, int n);
double norm1(double *x, int n);
double norm_diff (double *x, double *y, int n);
double norm1_diff(double *x, double *y, int n);

int mtxdet(double **c, int p, double *det);
void mtxswp(double **a, int p, int k);
void mtxinv(double **a, int p, double *det);
void mtxtra(double **a, double **b, int n, int p);
void vectra(double *a, double *b, int n);
void mtxmsc(double **a, double s, int n, int p);
double **mtxalloc(int n, int p);
void mtxfree(double **a, int n, int p);

void resdis(double **x, int n, int p, double *mu, double **sigma, double *d);
double rhobw(double u, double cc);
double mean_rhobw(double *u, double scale, int n, double cc);
double lossS(double *u, double scale, int n, double cc);
double scaleS(double *u, double kp, double cc, double initial_sc, int n);
void scaledpsi(double *u, double scale, int n, double cc, double *w);
int refine_s(double **x, int n, int p, double *init_mu, double **init_sigma, double init_scale, 
	     int k, int conv, double kp, double cc, 
	     double *mu, double **sigma, double *scale, double *rdis, double *weights);

/*
 * Compute FAST-S estimates of location and covariance matrix - 
 * similar to the one proposed for regression by 
 * Salibian-Barrera, M. and Yohai, V.J. (2005), 
 * "A fast algorithm for S-regression estimates". 
 * This current C implemention is by Valentin Todorov.
 *
 * INPUT:
 * X 	   - input data matrix - linear (FORTRAN) nn by pp
 * nsamp   - number of initial samples to generate
 * cc1     - constant for the biweight function
 * kkp     - constant (b) for the S-equation
 * bbest_r - number of 'best' solution to keep and iterate to convergence later
 * kstep
 * maxit
 * rtol
 *
 * OUTPUT:
 * ctr     - location vector
 * ccov    - covariance matrix (linear, p by p)
 * scale   - scale
 * converged
 * 
 */ 
void sest(double *X,
	  int *nn, int *pp, int *nsamp,
	  double *ctr, double *ccov, double *scale,
	  double *cc1, double *kkp,
	  int *bbest_r, int *kstep, int *maxit, double *rtol,
	  int *converged);

/* for "tracing" only : */
void disp_mat(double **a, int n, int m);
void disp_dble(double *a, int n);
void disp_int(int *a, int n);
void disp_lmat(double *a, int n, int p);

#define _use_blas_
/* ||x||_2 */
double norm(double *x, int n)
{
#ifdef  _use_blas_
    int one = 1;
    return(F77_CALL(dnrm2)(&n, x, &one));
#else
    int i;
    double s = 0.;
    for(i=0; i < n; i++)
	s += x[i] * x[i];
    return(sqrt(s));
#endif
}

double norm1(double *x, int n)
{
#ifdef  _USE_BLAS_
    int one = 1;
    return(F77_CALL(dasum)(&n, x, &one));
#else
    int i; double s = 0.;
    for(i=0; i < n; i++)
	s += fabs(x[i]);
    return(s);
#endif
}

/* ||x-y||_2 */
double norm_diff(double *x, double *y, int n)
{
    int i;
    double s = 0;
    for(i=0; i < n; i++)
	s += (x[i]-y[i])*(x[i]-y[i]);
    return( sqrt(s) );
}

/* ||x-y||_1 */
double norm1_diff(double *x, double *y, int n)
{
    int i;
    double s = 0;
    for(i=0; i < n; i++)
	s += fabs(x[i]-y[i]);
    return(s);
}

/* 
 * Resample without replacement:
 * Interface to R - call it from R instead of sample() for tests
 * when you want to have exactly the same sequence.
 *
*/
void r_sample(int *x, int *n, int *k)
{
    int i;
    resample(x, *n, *k);
    for(i=0; i<*k; i++)
	x[i] += 1;    
}

/* 
 * Resample without replacement - k out of n elements of the array 'array'
 *
*/
void resample(int *array, int n, int k)
{
   int i;
   for(i=0; i<n; ++i)  
       array[i] = i;
   selectwr(array, n, k);
   reverse(array, n);
}

void reverse(int *a, int n)
{
    int i, j, tmp;
    for(i=0, j=n-1; i<j; i++, j--)
    {
	tmp = a[i]; a[i] = a[j]; a[j] = tmp;
    }
}

/* Select from array without replacement */
void selectwr(int *array, int size, int nchoices)
{
   int i, temp, index;
   for(i=0; i<nchoices; ++i)
   {
       index = (double)size * unif_rand();
       temp = array[--size];
       array[size] = array[index];
       array[index] = temp;
   }
}

/*
 * best_r -  number of best solutions (out of nsamp) to keep
 * kstep  -  number of refining steps for each candidate
 * maxit  -
 * rtol   - 
 * converged - 
 */ 
void sest(double *X,
	  int *nn, int *pp, int *nsamp,
	  double *ctr, double *ccov, double *sscale,
	  double *cc1, double *kkp,
	  int *bbest_r, int *kstep, int *maxit, double *rtol,
	  int *converged)
{
    int i, j, ntry, ires, ind, iter, ibest;
    Rboolean singular, conv = 0;
    int n = *nn, p = *pp, nind = *pp+1, rank;
    int best_r = *bbest_r;
    double kp = *kkp, c1 = *cc1;
    double det;
    double scale, s_test, s_best, s_worst=INFI, s_tmp;

    /* (Pointers to) Arrays - to be allocated */
    int *b_i   = (int *) R_alloc(n, sizeof(int));
    double *mu = (double *) R_alloc(p, sizeof(double));
    double *mu_best = (double *) R_alloc(p, sizeof(double));
    double **x 	      = mtxalloc(n, p);
    double **tmp_mat  = mtxalloc(p, p);
    double **sigma    = mtxalloc(p,p);
    double **sigma_best = mtxalloc(p,p);
    double *rdis    = (double *) R_alloc(n, sizeof(double));
    double *weights = (double *) R_alloc(n, sizeof(double));
    double *tmp     = (double *) R_alloc(n, sizeof(double));		
   
    double *best_scales = (double *) R_alloc(best_r, sizeof(double));		
    double *best_mus = (double *) R_alloc(best_r*p, sizeof(double));		
    double *best_sigmas = (double *) R_alloc(best_r*p*p, sizeof(double));		
    for(i=0; i < best_r; i++) 
	best_scales[i] = INFI;

    // Copy the linear (FORTRAN) data matrix into a C matrix[][].
    for(i=0; i < n; i++)
	for(j=0; j < p; j++)
	    x[i][j] = X[j*n+i];

//    disp_mat(x, n, p);
//    disp_lmat(X, n, p);
//    disp_dble(X, n*p);
//    disp_int(&best_r,1);
    
    GetRNGstate();

    for(ires=0; ires < *nsamp; ires++) 
    {

	/* find a candidate 
	 *
	 * Take a sample of p+1 indices and compute the mean,
	 * covariance matrix and determinant of the resulting
	 * subsample. If the matrix is singular, retry until
	 * max number of retrieals reached.
	 */
	ntry = 0;		// maximum number of singular p+1
				// subsamples before giving up 
	singular = TRUE;
	while(singular) 
	{
	    R_CheckUserInterrupt();
	    if(++ntry > MAX_NTRY) 
	    {
		REprintf("\nToo many singular resamples\nAborting sest()\n\n");
		*sscale = -1.0;
		goto cleanup_and_return;
	    }

	    resample(b_i, n, p+1);
	    //covp(x, n, p, b_i, p+1, tmp, tmp_mat);
	    //singular = mtxdet(tmp_mat, p, &det);
            covp(x, &n, &p, b_i, &nind, tmp, tmp_mat, &det, &rank); 
	    singular = rank < p;
//  	    REprintf("\n%d: Singular=%d, LOGDET= %lf\n", ires, singular, det);
	}

	mtxmsc(tmp_mat, pow(exp(det), -1.0/p), p, p);

        // Perform k steps of iterative reweighting on the elemental set
        if(*kstep > 0)
        {
            // do the refining: no convergence, only several (kstep) iterations
	    refine_s(x, n, p, tmp, tmp_mat, -1, *kstep, conv, kp, c1,
			mu, sigma, &scale, rdis, weights);
        } else
        {
            // k = 0 means "no refining"
            vectra(tmp, mu, p);
            mtxtra(tmp_mat, sigma, p, p);
            resdis(x, n, p, mu, sigma, rdis);
            scale = mymedabs(n, rdis)/.6745;
        }

        if(ires > 0)
        {
            /* 
	     * if this isn't the first iteration....
             * check whether new mu/sigma belong to the top best results; if so keep
             * mu and sigma with corresponding scale.
	     */
	    s_test = lossS(rdis, s_worst, n, c1);
	    if(s_test < kp)
	    {
	        s_best = scaleS(rdis, kp, c1, scale, n);
		ind = maxind(best_scales, &s_worst, best_r);
		best_scales[ind] = s_best;
	    	for(i=0; i<p; i++)
	    	{
		    best_mus[ind*p+i] = mu[i];
	            for(j=0; j<p; j++)
		        best_sigmas[ind*p*p+i*p+j] = sigma[i][j];
	    	}
		maxind(best_scales, &s_worst, best_r);
  	    	
//		REprintf("\n%d: Best_r=%d, Scale= %lf20.16\n", ires, best_r, s_test);
//		disp_mat(sigma, p, p);

	    }
        } else
        {
            /*
	     * if this is the first iteration, then this 
	     * is the best solution anyway... store it.
	     */
	    ind = best_r-1;
	    best_scales[ind] = scaleS(rdis, kp, c1, scale, n);
	    for(i=0; i<p; i++)
	    {
		best_mus[ind*p+i] = mu[i];
	        for(j=0; j<p; j++)
		    best_sigmas[ind*p*p+i*p+j] = sigma[i][j];
	    }
        }
    }
    
//    Rprintf("\nNow refine() to convergence for %d very best ones:\n", best_r);
//    disp_dble(best_scales, best_r);
//    disp_lmat(best_mus, best_r, p);
//    disp_lmat(best_sigmas, best_r*p, p);


    /* Now look for the very best:
     * Do the complete refining step until convergence (conv=1) starting
     * from the best subsampling candidates (possibly refined)
     */
    s_best = INFI; 
    for(ibest=0; ibest<best_r; ibest++) 
    {
	conv = 1;

	// Extract the ind-th set if mu, sigma and scale
	s_tmp = best_scales[ibest];
	for(i=0; i<p; i++)
	{
	    tmp[i] = best_mus[ibest*p+i];
            for(j=0; j<p; j++)
   	        tmp_mat[i][j] = best_sigmas[ibest*p*p+i*p+j];
	}
        
	iter = refine_s(x, n, p, tmp, tmp_mat, s_tmp, 1, conv, kp, c1,
			mu, sigma, &scale, rdis, weights);

	if(scale < s_best) 
	{
	    // we are improving - remember the best so far estimates
	    s_best = scale;
	    vectra(mu, mu_best, p);
	    mtxtra(sigma, sigma_best, p, p);
	}

	if(!conv && *converged) 
	    *converged = 0;
	if(*kstep < iter) 
	    *kstep = iter;
    }

    mtxmsc(sigma_best, s_best*s_best, p, p);
    
    *maxit = *kstep;
    *sscale = s_best;
    vectra(mu_best, ctr, p);
    for(i=0; i<p; i++)
	for(j=0; j<p; j++)
	    ccov[j*p+i] = sigma_best[i][j];
    
cleanup_and_return:

    PutRNGstate();

    mtxfree(x, n, p);
    mtxfree(tmp_mat, p, p);
    mtxfree(sigma, p, p);
    mtxfree(sigma_best, p, p);
}

/*
 * 	Do "k" iterative reweighting refining steps from 
 * 	initial mu and sigma
 *
 * 	init_scale - if present (>= 0), it's used, o/w the MAD is used
 *	k = number of refining steps
 *   	conv = 0 means "do k steps and don't check for convergence"
 *   	conv = 1 means "stop when convergence is detected, or the
 *   		maximum number of iterations is achieved"
 *   	kp and cc = tuning constants of the equation
 */
int refine_s(double **x, int n, int p, double *init_mu, double **init_sigma, double init_scale, 
	     int k, int conv, double kp, double cc, 
	     double *mu, double **sigma, double *scale, double *rdis, double *weights)
{
    double convtol = 1e-20;
    double *mu_1 = (double *)R_Calloc(p, double);
    double **sigma_1 = mtxalloc(p,p);
    double **tmp_mat = mtxalloc(p,p);
    double sc, det;
    int i, singular;

    resdis(x, n, p, init_mu, init_sigma, rdis);
    if(init_scale < 0.0)
        init_scale = sc = mymedabs(n, rdis)/.6745;
    else
	sc = init_scale;

    // if conv == 1 then set the max no. of iterations to 50 
    // magic number alert!!!
    if(conv == 1)
        k = 50;

    // copy the initial estimates mu and sigma to work storage
    vectra(init_mu, mu, p);
    mtxtra(init_sigma, sigma, p, p);
    for(i=0; i<k; i++)
    {
        // do one step of the iterations to solve for the scale
        sc = sqrt(sc*sc * mean_rhobw(rdis, sc, n, cc) / kp);

        // now do one step of reweighted mean and covariance 
	// with the "improved scale"
 	scaledpsi(rdis, sc, n, cc, weights);
	covwt(x, n, p, weights, mu_1, sigma_1);
	singular = mtxdet(sigma_1, p, &det);
//  	REprintf("\n------Singular=%d, LOGDET= %lf\n", singular, det);

	if(singular)
	{
	    vectra(init_mu, mu_1, p);
	    mtxtra(init_sigma, sigma_1, p, p);
	    *scale = init_scale;
	    break;
	}

	mtxmsc(sigma_1, pow(exp(det),-1.0/p), p, p);
        if(norm_diff(mu, mu_1, p) / norm(mu, p) < convtol)
            break;

	vectra(mu_1, mu, p);
	mtxtra(sigma_1, sigma, p, p);
	resdis(x, n, p, mu_1, sigma_1, rdis);
    
//	disp_dble(mu, p);
//    	disp_mat(sigma, p, p);
//    	disp_dble(&sc, 1);
    }
    *scale = sc;
	
    R_Free(mu_1);
    mtxfree(sigma_1,p,p);
    mtxfree(tmp_mat,p,p);

    return (i);	
}

/*
 * Compute the location vector and covariance matrix of a subsample
 * of the data matrix a[n][p] with n rows and p columns. The subsample is 
 * given by the vector id with length p.
 */ 
void covpold(double **a, int n, int p, int *id, int np, double *t, double ** cov)
{
    int i, j, k;
    for(j=0; j<p; j++)
    {
        t[j] = 0.0;
	for(k=0; k<p; k++)
	    cov[j][k] = 0.0;
    }
    for(i=0; i<np; i++)
    {
	for(j=0; j<p; j++)
	{
	    t[j] += a[id[i]][j];
	    for(k=0; k<=j; k++)
		cov[j][k] += a[id[i]][j]*a[id[i]][k];
	}
    }
    for(j=0; j<p; j++)
    {
	for(k=0; k<=j; k++)
	{
	    cov[j][k] -= t[j]*t[k]/np;
	    cov[j][k] /= np-1;
	    cov[k][j]  = cov[j][k];
	}
    }

    for(j=0; j<p; j++)
	t[j] /= np;
}

/*
 * Compute the location vector and covariance matrix of a subsample
 * of the data matrix a[n][p] with n rows and p columns.
 * The covariance is computed using QR decomposition. The subsample is
 * given by the vector indices with length nind. Return the mean and cov[][].
 * Compute the rank and if full rank, compute the determinant.
 */
void covp(double **x, int *n, int *p, int *indices, int *nind, 
		double *mean, double **cov, double *det, int *rank)
{
int i, j, k;
int pp = *p, nnind = *nind;
double tol = 1e-7, s=0;

double *cx = (double *) R_alloc(pp*pp, sizeof(double));		
double *xw = (double *) R_alloc(nnind*pp, sizeof(double));		
double *qraux = (double *) R_alloc(pp, sizeof(double));		
double *work = (double *) R_alloc(2*pp, sizeof(double));		
int *pivot = (int *) R_alloc(pp, sizeof(int));

    // compute the mean, put submatrix into xw, center its columns
    for(j=0; j<pp; j++) 
    {
        mean[j] = 0.0;
        for(i=0; i<nnind; i++) {
            mean[j] += (xw[i+j*nnind] = x[indices[i]][j]) / (double) nnind;
        }
        for(i=0; i<nnind; i++) {
            xw[i+j*nnind] -= mean[j];
        }
    }

   // QR decomposition of the submatrix
   F77_CALL(dqrdc2)(xw, nind, nind, p, &tol, rank, qraux, pivot, work);

   // build the cov matrix of the subsample using the QR decomp
   for(i=0; i<pp; i++)
      for(j=i; j<pp; j++) 
      {
         s = 0;
	 for(k=0; k<=i; k++) 
	    s += xw[k + j*nnind] * xw[k + i*nnind];
	 cx[j + i*pp] = cx[i + j*pp] = (s / (double) (nnind - 1));
      }

   for(i=0; i<pp; i++)
      for(j = 0; j < pp; j++)
      cov[i][j] = cx[i + pp*j];

   /* if full rank, compute det of cov matrix */
   /* det^2 = (nind-1)^p * det(cov matrix) */
   if(*rank == pp) 
   { 
      *det = 1.0;
      for(j=0; j<pp; j++)
         *det += log(fabs(xw[j + nnind*j]));
   } else
	*det = log(0); 
}



/*
 * Compute the location vector and covariance matrix of a 
 * data matrix a[n][p] with n rows and p columns.
 */ 
void covar(double **a, int n, int p, double *t, double ** cov)
{
    int i;	
    double det;
    int rank;
    int *id = (int *)R_Calloc(n, int);
    for(i=0; i<n; i++)
	id[i] = i;
    //covp(a, n, p, id, n, t, cov);
    covp(a, &n, &p, id, &n, t, cov, &det, &rank); 
    R_Free(id);
}

/*
 * Compute weighted mean and unscaled weighted 
 * covariance matrix from the
 * observations a[][] and using weights in wts[].
 */ 
void covwt(double **a, int n, int p, double *wts, double *t, double **cov)
{
    int i, j, l;
    double sumw=0, sc;
    for(i=0; i<n; i++)
	sumw += wts[i];
    for(j=0; j<p; j++)
    {
	sc = 0;
	for(i=0; i<n; i++)
	    sc += a[i][j]*wts[i];
	t[j] = sc/sumw;
    }
    
    for(i=0; i<p; i++)
    {
	for(j=0; j<=i; j++)
	{
	    sc = 0.0;
	    for(l=0; l<n; l++)
                sc += wts[l]*(a[l][i] - t[i])*(a[l][j] - t[j]);
	    cov[i][j] = sc;
	    cov[j][i] = cov[i][j];
	}
    }
}

/*
 * Returns square root of the mahalanobis distances of x 
 * with respect to mu and sigma.
 */
void resdis(double **x, int n, int p, double *mu, double **sigma, double *d)
{
    int l, i, j;
    double temp;
    double **sigma_1 = mtxalloc(p,p);
    
    for(i=0; i<n; i++)
        d[i] = 0.0;

    mtxtra(sigma, sigma_1, p, p);
    mtxinv(sigma_1, p, &temp);
    if(temp >= 0)	// not singular
    {
        for(l=0; l<n; l++)
        {
            d[l] = 0.0;
	    for(i=0; i<p; i++)
	    {
	        temp = x[l][i] - mu[i];
	        for(j=0; j<p; j++)
		    d[l] += temp * (x[l][j] - mu[j]) * sigma_1[i][j];
	    }
	    d[l] = sqrt(d[l]);
        }
    }

    mtxfree(sigma_1, p, p);
}

void mtxswp(double **a, int p, int k)
{
    double b, d;
    int i, j;

    d = a[k][k];
    for(j=0; j< p; j++)
        a[j][k] /= d;
    for(i=0; i<p; i++)
    {
	if(i != k)
	{
	    b = a[k][i];
	    for(j=0; j<p; j++)
	        a[j][i] -= b*a[j][k];
	    a[k][i] = -b/d;
	}
    }
    a[k][k] = 1/d;
}

void mtxinv(double **a, int p, double *det)
{
    double pivot;
    int j;
    double eps;

    if(p < 5)               eps = 1.0E-12;
    else if(p > 5 && p < 8) eps = 1.0E-14;
    else            	    eps = 1.0E-16;

    *det = 1.0;
    for(j=0; j<p; j++)
    {
	pivot = a[j][j];
	*det *= pivot;
	if(pivot < eps)
	{
	    det = 0;
	    return;
	}
	mtxswp(a, p, j);
    }
}

// b <- a
void mtxtra(double **a, double **b, int n, int p)
{
    int i, j;

    for(i=0; i<n; i++)
	for(j=0; j<p; j++)
	    b[i][j] = a[i][j];
}

// Allocate storage for an nxp matrix
double **mtxalloc(int n, int p)
{
    int i;

    double **a  = (double **)R_Calloc(n, double *);		
    for(i=0; i<n; i++) 
	a[i] = (double *)R_Calloc(p, double);
    return a;
}

// Free the allocated  storage for an nxp matrix
void mtxfree(double **a, int n, int p)
{
    int i;
    for(i=0; i<n; i++) 
	R_Free(a[i]);
    R_Free(a);
}

// b <- a
void vectra(double *a, double *b, int n)
{
    int i;

    for(i=0; i<n; i++)
        b[i] = a[i];
}

// a <- s * a
void mtxmsc(double **a, double s, int n, int p)
{
    int i, j;
    for(i=0; i<n; i++)
	for(j=0; j<p; j++)
	    a[i][j] *= s;
}

double mymed(int n, double *x)
{
  int k = -1;

  if(n%2) {
    k = (n-1)/2;
    rPsort(x, n, k);
    return(x[k]);
  }

  k = n/2;
  rPsort(x, n, k);
  rPsort(x, k, k-1);
  return((x[k-1] + x[k]) / 2.0);
}

double mymedabs(int n, double *x)
{
    int i;
    double tmp;
    double *vec = (double *)R_Calloc(n, double);
    for(i=0; i<n; i++)
	vec[i] = fabs(x[i]);

    tmp = mymed(n, vec);
    R_Free(vec);
    return(tmp);
}

double mad(int n, double *x, double *dwork1, double *dwork2, double *mu)
{
  const int IONE = 1;
  int i = 0;
  F77_CALL(dcopy)(&n, x, &IONE, dwork1, &IONE);
  *mu = mymed(n, dwork1);
  for(i = 0; i < n; i++)
    dwork1[i] = fabs(dwork1[i] - *mu);
  return mymed(n, dwork1) * 1.4826;
}

/*
 * Find the maximal value inthe arra a[] and return the 
 * corresponding index.
 */
int maxind(double *a, double *maxval, int n)
{
    if(n == 1)
    {
    	*maxval = a[0];
	return(0);
    } else
    {
	int i, k = 0;
	double tt = a[0];
	for(i=1; i < n; i++)
	{
	    if(tt < a[i]) 
	    {
		tt = a[i];
		k = i;
	    }
	}
	*maxval = tt;
	return(k);
    }
}


/*
The bisquare rho function:

              |   x^2/2 - x^4/2*c1^2 + x^6/6*c1^4         |x| <=  c1
  rho(x) =    |
              |   c1^2/6                                  |x| > c1

*/
double rhobw(double u, double cc)
{
    double t2, c2 = cc*cc;
    double ret = c2/6;

    if(fabs(u) <= cc)
    {
	t2 = u*u;
	ret = t2*(1.0 - t2/c2 + t2*t2/(3*c2*c2))/2;
    }
    return ret;
}
	
/*
 * Return the mean of the rho(u/scale) for the Tukey's
 * biweight function rho()
 *	1/n * sum(rho(u/s))
 */
double mean_rhobw(double *u, double scale, int n, double cc)
{
    int i;
    double s = 0;
    for(i=0; i < n; i++)
	s += rhobw(u[i]/scale, cc);
    return (s/n);
}


/*
 * The objective function:
 * we solve loss.S(u, s, cc) = b for "s"
 */
double lossS(double *u, double scale, int n, double cc)
{
    return mean_rhobw(u, scale, n, cc);
}

/*
 * Scaled psi function for Tukey's biweight to comute the weights:
 * Weight function = psi(u)/u
 */ 
void scaledpsi(double *u, double scale, int n, double cc, double *w)
{
    int i=0;
    double t;
    // ##ifelse(abs(xx) < c1, xx - 2 * xx^3/c1^2 + xx^5/c1^4, 0)
    for(i=0; i<n; i++)
    {
	w[i] = 0;
	t = u[i]/scale;
//
//	VT::25.08.2016 - warnung reported by Duncan Murdoch
//		for 1.4.1
//		
//	if(fabs(t <= cc))
	if(fabs(t) <= cc)
	{
	    t = 1 - t*t/cc/cc;
	    w[i] = t * t;
	    w[i] *= cc*cc/6.0;
	}
    }
}

/*
 * Computes Tukey's biweight objective function (scale)
 * (respective to the mahalanobis distances u) using the
 *  rho() function and the konstants kp and c1
 */
double scaleS(double *u, double kp, double cc, double initial_sc, int n)
{
    // find the scale, full iterations
    int maxit = 200, it = 0;
    double sc = initial_sc, sc2;
    double eps = 1e-20;
    double err = 1.0;

    while(++it < maxit && err > eps)
    {
	sc2 = sqrt(sc*sc * mean_rhobw(u, sc, n, cc) / kp);
	err = fabs(sc2/sc - 1.0);
	sc = sc2;
    }

    return(sc);
}

void disp_dble(double *a, int n)
{
    int i;
    for(i=0; i < n; i++) Rprintf("%lf ",a[i]);
    Rprintf("\n");
}

void disp_int(int *a, int n)
{
    int i;
    for(i=0; i < n; i++) 
	Rprintf("%d ", a[i]);
    Rprintf("\n");
}

void disp_mat(double **a, int n, int m)
{
    int i,j;
    for(i=0; i < n; i++) {
	Rprintf("\n");
	for(j=0; j < m; j++) Rprintf("%10.8f ",a[i][j]);
    }
    Rprintf("\n");
}

/*
 * Display an nxp matrix stored in linear array in 
 * row major order (C)
 */ 
void disp_lmat(double *a, int n, int p)
{
    int i,j;
    for(i=0; i < n; i++) 
    {
	Rprintf("\n");
	for(j=0; j < p; j++) 
	    Rprintf("%10.8f ", a[i*p+j]);
    }
    Rprintf("\n");
}

/*
   Compute the determinant and check if the matrix 'c' is singular
   using QR decomposition
*/
int mtxdet(double **c, int p, double *det)
{

    int i, j, k;
    int rank;
    double sum, tol = 1.0e-7;

    double *cx = (double *) R_alloc(p*p, sizeof(double));		
    double *qraux = (double *) R_alloc(p, sizeof(double));		
    double *work = (double *) R_alloc(2*p, sizeof(double));		

    int *pivot = (int *) R_alloc(p, sizeof(int));

    for(i=0; i<p; i++)
        for(j = 0; j < p; j++)
	    cx[i + p*j] = c[i][j];

    F77_CALL(dqrdc2)(cx, &p, &p, &p, &tol, &rank, qraux, pivot, work);
    if(rank < p) 
        return (1);

    sum = 0.0;
    for(k = 0; k < p; k++)
	sum += log(fabs(cx[k + p*k]));
    *det = sum;

    return (0);
}
 

