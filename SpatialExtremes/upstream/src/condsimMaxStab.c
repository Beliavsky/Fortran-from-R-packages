#include "header.h"

void stirling2ndKind(int *n, int *k, double *ans){
  /* This function computes the Stirling numbers of the second kind,
     i.e., the number of partition of a set with n elements of size
     k */

  double a[*n + 1][*k + 1];

  for (int i=0;i<=*n;i++)
    for (int j=0;j<=*k;j++)
      a[i][j] = 0;

  for (int i=0; i<=*n; i++){
    for (int j=0;j<=*k;j++){
      if (i < j)
	a[i][j] = 0;

      else if (i == j)
	a[i][j] = 1;

      else if (j == 1)
	a[i][j] = 1;

      else if (j == 0)
	a[i][j] = 0;

      else if (j<i)
	a[i][j] = j * a[i-1][j] + a[i-1][j-1];

    }
  }

  *ans = a[*n][*k];
  return;
}

void bell(int *n, int *ans){

  /* This function computes the Bell numbers, i.e., the number of
     partition of a set with n elements. */

  int bell[*n + 1];
  bell[0] = 1;

  int iter = 1;
  while (iter <= *n){
    bell[iter] = 0;

    for (int i = 0; i < iter; i++)
      bell[iter] += choose(iter-1, i) * bell[i];

    iter++;
  }

  *ans = bell[*n];
  return;
}

void listAllPartOfASet(int *n, int *nPart, int *allPart, int *allSize){
  /* This function lists all the partitions of a set with n
     elements. It is based on a paper of Michael Orlov [2002].

  The convention used here is the following. The partition {{x_1,
  x_3}, {x_2, x_5}, {x_4, x_6}} is identified to (0, 1, 0, 2, 1,
  2). */

  int k[*n], m[*n];

  // Initialize k
  for (int i=0;i<*n;i++)
    allPart[i] = k[i] = 0;

  allSize[0] = 1;

  unsigned int currentPart = 0;

  while(currentPart < *nPart){
    // Restore m as it is invariant
    m[0] = k[0];
    for (unsigned i=1;i<*n;i++)
      m[i] = imax2(m[i-1], k[i]);

    currentPart++;

    for (unsigned i = (*n - 1); i > 0; --i){
      if (k[i] <= m[i-1]){
	++k[i];

	const int newMax = imax2(m[i], k[i]);
	m[i] = newMax;

	for (unsigned j = (i+1); j < *n; ++j){
	  k[j] = 0;
	  m[j] = newMax;
	}

	for (unsigned j=0;j<*n;j++)
	  allPart[currentPart * *n + j] = k[j];

	allSize[currentPart] = m[*n - 1] - m[0] + 1;

	break;
      }
    }
  }

  return;
}

void pmvt(double *bounds, int *n, double *DoF, double *mu, double *scaleMat,
	  double *prob, double *err, int *nMc){
  /* This function computes the CDF of a multivariate Student
     distribution random vector using quasi monte carlo sequence*/

  const int primeNumbers[100] = {
    2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
    31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
    73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
    353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
    419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
    467, 479, 487, 491, 499, 503, 509, 521, 523, 541};

  // 0. Centering the data
  double *boundsc = malloc(*n * sizeof(double));
  for (int i=0;i<*n;i++)
    boundsc[i] = bounds[i] - mu[i];

  // 1. We need to reorder the variables
  double *newCor = malloc(*n * *n * sizeof(double));
  memcpy(newCor, scaleMat, *n * *n * sizeof(double));

  //Compute the Cholesky decomposition of newCor
  int info = 0;
  F77_CALL(dpotrf)("L", n, newCor, n, &info FCONE);

  if (info != 0)
    error("1. error code %d from Lapack routine '%s'", info, "dpotrf");

  // Get the inverse entries of newCor for efficiency
  double *inewCor = malloc(*n * *n * sizeof(double));
  for (int i=0;i<*n * *n;i++)
    inewCor[i] = 1 / newCor[i];

  const double isqrtDoF = 1 / sqrt(*DoF);

  // Initialize objects
  int M = 1, N = *n * 1000;
  double T = 0, E = 1, V = 0;
  double *I = malloc(M * sizeof(double)),
    *Delta = malloc(*n * sizeof(double)),
    delta,
    *q = malloc(*n * sizeof(double)),
    *e = malloc(*n * sizeof(double)),
    *f = malloc(*n * sizeof(double)),
    *w = malloc(*n * sizeof(double)),
    *y = malloc((*n - 1) * sizeof(double)),
    *ea = malloc(*n * sizeof(double)),
    *fa = malloc(*n * sizeof(double)),
    *wa = malloc(*n * sizeof(double)),
    *ya = malloc((*n - 1) * sizeof(double)),
    *b = malloc(*n * sizeof(double)),
    *ba = malloc(*n * sizeof(double));//a for antithethic

  for (int i=0;i<*n;i++)
    q[i] = sqrt(primeNumbers[i]);

  GetRNGstate();

  for (int i=0;i<M;i++){
    I[i] = 0;

    for (int j=0;j<*n;j++)
      Delta[j] = unif_rand();

    for (int j=0;j<N;j++){
      for (int k=0;k<*n;k++){
	double dummy = j * q[k] + Delta[k];
	w[k] = fabs(2 * (dummy - ftrunc(dummy)) - 1);
	wa[k] = 1 - w[k];
      }

      double s = sqrt(qchisq(w[*n-1], *DoF, 1, 0)),
	sa = sqrt(qchisq(wa[*n-1], *DoF, 1, 0));

      for (int k=0;k<*n;k++){
	b[k] = boundsc[k] * s * isqrtDoF;
	ba[k] = boundsc[k] * sa * isqrtDoF;
      }

      e[0] = pnorm(b[0] * inewCor[0], 0, 1, 1, 0);
      ea[0] = pnorm(ba[0] * inewCor[0], 0, 1, 1, 0);
      f[0] = e[0];
      fa[0] = ea[0];

      for (int k=1;k<*n;k++){
	y[k-1] = qnorm(w[k-1] * e[k-1], 0, 1, 1, 0);
	ya[k-1] = qnorm(wa[k-1] * ea[k-1], 0, 1, 1, 0);

	double dummy = 0, dummya = 0;
	for (int l=0;l<k;l++){
	  dummy += newCor[k + l * *n] * y[l];
	  dummya += newCor[k + l * *n] * ya[l];
	}

	e[k] = pnorm((b[k] - dummy) * inewCor[k * (*n + 1)], 0, 1, 1, 0);
	ea[k] = pnorm((ba[k] - dummya) * inewCor[k * (*n + 1)], 0, 1, 1, 0);
	f[k] = e[k] * f[k - 1];
	fa[k] = ea[k] * fa[k - 1];
      }

      I[i] = I[i] + (0.5 * (f[*n - 1] + fa[*n - 1]) - I[i]) / ((double) j + 1);
    }

    delta = (I[i] - T) / ((double) i + 1);
    T += delta;
    V = (i - 1) * V / ((double) i + 1) + delta * delta;
    E = 3 * sqrt(V);
  }

  PutRNGstate();

  *prob = T;
  *err = E;
  *nMc = (int )N;

  free(newCor); free(inewCor); free(I); free(Delta);
  free(q); free(e); free(f); free(w); free(y);
  free(ea); free(fa); free(wa); free(ya);
  free(b);
  free(ba);
  free(boundsc);
  return;
}

void computeWeightsSC(int *nCond, double *y, int *nPart, int *allPart,
		    int *allSize, double *cov, double *weights){
  /* This function computes the distribution function of the random
     partition given the conditioning values, i.e., Pr[\theta = \tau |
     Z(x) = z]. (Schlather case) */

  int *currentPart = malloc(*nCond * sizeof(int));

  // Set the weights to one
  for (int i=0;i<*nPart;i++)
    weights[i] = 1;

  for (int i=0;i<*nPart;i++){
    // Loop over all the partitions
    memcpy(currentPart, allPart + i * *nCond, *nCond * sizeof(int));

    for (int j=0;j<allSize[i];j++){
      //Loop over all sets of the current partition

      int r = 0;
      for (int k=0;k<*nCond;k++)
	r += (currentPart[k] == j);

      int nr = *nCond - r,
	*tau = malloc(r * sizeof(int)),
	*taubar = malloc(nr * sizeof(int));

      gettau(nCond, currentPart, &j, tau);

      if (r < *nCond) {
	gettaubar(nCond, currentPart, &j, taubar);

	/* Define the parameters of the distribution of the extremal
	   function. For the Schlather case, it is a multivariate
	   Student distribution. Hence below we compute the degree of
	   freedom (DoF), the Cholesky decomposition of the scale
	   matrix (scaleMat) and the location parameter (mu) */
	double DoF = 0,
	  *scaleMat = malloc(nr * nr * sizeof(double)),
	  *mu = malloc(nr * sizeof(double));
	getParametersSC(tau, taubar, &r, &nr, cov, y, &DoF, mu, scaleMat);

	/* Compute the proba values, i.e., the probabilities that one
	   extremal function is smaller than the hitting bounds at the
	   remaining conditioning locations. */
	double prob = 0;
	computeprobaSC(&DoF, mu, scaleMat, y, &nr, taubar, &prob);
	weights[i] *= prob;

	free(scaleMat); free(mu);
      }

      /* Compute the intensity function values, i.e.,
	 \lambda_{x_{\tau_j}}(z_{\tau_j}). */
      double f = 0;
      getfvaluesSC(y, nCond, &r, tau, cov, &f);

      weights[i] *= exp(f);
      free(tau); free(taubar);
    }
  }

  /* Compute the normalizing constant and normalize the
     distribution */
  double isumWeights = 0;
  for (int i=0;i<*nPart;i++)
    isumWeights += weights[i];

  isumWeights = 1 / isumWeights;

  for (int i=0;i<*nPart;i++)
    weights[i] *= isumWeights;

  free(currentPart);
  return;
}

void computeprobaSC(double *DoF, double *mu, double *scaleMat, double *y,
		  int *ntaubar, int *taubar, double *prob){

  *prob = 0;

  double *ytaubar = malloc(*ntaubar * sizeof(double));
  for (int i=0;i<*ntaubar;i++)
    ytaubar[i] = y[taubar[i]];

  if (*ntaubar == 1)
    *prob = pt(*ytaubar - *mu, *DoF, 1, 0);

  else {
    /* Computing the proba */
    double err = 0;
    int nMc = 0;

    pmvt(ytaubar, ntaubar, DoF, mu, scaleMat, prob, &err, &nMc);
  }

  free(ytaubar);

  return;
}

void getfvaluesSC(double *y, int *n, int *ntau, int *tau, double *cov,
		double *f){

  /* This function computes the logarithm of the intensity
     function, i.e., \lambda_{x_\tau_j}(z_{\tau_j}). */

  /* Get the sub-matrix of cov */
  double *covjchol = malloc(*ntau * *ntau * sizeof(double));
  getSubMatrix(cov, n, ntau, tau, ntau, tau, covjchol);

  /* Get the Cholesky decomposition of this submatrix */
  int info = 0;
  F77_CALL(dpotrf)("U", ntau, covjchol, ntau, &info FCONE);

  if (info != 0)
    error("4. error code %d from Lapack routine '%s'", info, "dpotrf");

  //Get the sub-vector of y
  double *yj = malloc(*ntau * sizeof(double));
  for (int k=0;k<*ntau;k++)
    yj[k] = y[tau[k]];

  /* Compute the log of the determinant of cov */
  int oneInt = 1;
  double det = 0;
  for (int i=0;i<*ntau;i++)
    det += log(covjchol[i * (1 + *ntau)]);

  det *= 2;

  /* Compute covjchol^{-T} %*% y */
  double *covjCholMean = malloc(*ntau * sizeof(double));
  memcpy(covjCholMean, yj, *ntau * sizeof(double));

  F77_CALL(dtrsv)("U", "T", "N", ntau, covjchol, ntau, covjCholMean, &oneInt
       FCONE FCONE FCONE);

  double mahal = 0;
  for (int i=0;i<*ntau;i++)
    mahal += covjCholMean[i] * covjCholMean[i];

  *f = (1 - *ntau) * M_LN_SQRT_PI - 0.5 * det - 0.5 * (*ntau + 1) *
    log(mahal) + lgammafn(0.5 * (*ntau + 1));

  free(covjCholMean); free(covjchol);
  return;
}

void getSubMatrix(double *mat, int *dim, int *nr, int *rows, int *nc,
		  int *cols, double *submat){

  for (int i=0;i<*nr;i++)
    for (int j=0;j<*nc;j++)
      submat[i + *nr * j] = mat[rows[i] + *dim * cols[j]];

  return;
}

void getParametersSC(int *tau, int *taubar, int *ntau, int *ntaubar, double *cov,
		   double *y, double *DoF, double *mu, double *scaleMat){

  /* This function computes the parameters of distribution of the
     extremal functions. For the Schlather case, it is a multivariate
     Student distribution with degree of freedom (DoF), location
     parameter (mu) and scale matrix (scaleMat). */

  int oneInt = 1;

  // First get the submatrices of Sigma and the subvector of y
  int dim = *ntau + *ntaubar;

  double *covCholx = malloc(*ntau * *ntau * sizeof(double));
  getSubMatrix(cov, &dim, ntau, tau, ntau, tau, covCholx);

  int info = 0;
  F77_CALL(dpotrf)("U", ntau, covCholx, ntau, &info FCONE);

  if (info != 0)
    error("0. error code %d from Lapack routine '%s'", info, "dpotrf");

  double *covs = malloc(*ntaubar * *ntaubar * sizeof(double));
  getSubMatrix(cov, &dim, ntaubar, taubar, ntaubar, taubar, covs);

  double *covsx = malloc(*ntaubar * *ntau * sizeof(double));
  getSubMatrix(cov, &dim, ntaubar, taubar, ntau, tau, covsx);

  double *yx = malloc(*ntau * sizeof(double));
  for (int i=0;i<*ntau;i++)
    yx[i] = y[tau[i]];

  /* Get the degree of freedom parameter */
  *DoF = *ntau + 1;

  /* Compute the location parameter:

     It is given by mu = covsx %*% covx^(-1) %*% yx, so in terms of
     the Cholesky decomposition we will compute it as follows

     mu = covsx %*% (covCholx^T %*% covCholx)^(-1) %*% yx
        = covsx %*% covCholx^(-1) %*% covCholx^(-T) %*% yx,

     i.e.,
       a) solve the (triangular) system X %*% covCholx = covsx -->> dummy1
       b) solve the (triangular) system covCholx^T %*% Y = yx  -->> dummy2
       c) Then mu = X %*% Y.

     Note that Y is a vector and X is a matrix!

  */

  double *dummy1 = malloc(*ntaubar * *ntau * sizeof(double));
  memcpy(dummy1, covsx, *ntaubar * *ntau * sizeof(double));

  double alpha = 1;
  F77_CALL(dtrsm)("R", "U", "N", "N", ntaubar, ntau, &alpha, covCholx, ntau,
		  dummy1, ntaubar FCONE FCONE FCONE FCONE);

  double *dummy2 = malloc(*ntau * sizeof(double));
  memcpy(dummy2, yx, *ntau * sizeof(double));
  F77_CALL(dtrsv)("U", "T", "N", ntau, covCholx, ntau, dummy2, &oneInt
		  FCONE FCONE FCONE);

  double beta = 0;
  /*
    for (int i=0;i<*ntaubar;i++)
    mu[i] = 0;

    Not necessary according to the documentation of BLAS...

  */

  F77_CALL(dgemv)("N", ntaubar, ntau, &alpha, dummy1, ntaubar, dummy2, &oneInt,
		  &beta, mu, &oneInt FCONE);

  /* Compute the scale matrix:

     It is given by

     scaleMat = mahal / (k+1) * (covs - covsx %*% covx^(-1) %*% covsx^T),

     so in terms of the Cholesky decomposition we will compute it as
     follows

     dummy3 = covsx %*% covx^(-1) %*% covsx^T
            = covsx %*% (covCholx^T %*% covCholx)^(-1) %*% covsx^T
	    = covsx %*% covCholx^(-1) %*% covCholx^(-T) %*% covsx^T
	    = covsx %*% covCholx^(-1) %*% (covsx %*% covCholx^(-1))^T
	    = dummy1 %*% dummy1^T.

     dummy3 = covs - dummy3;

     mahal = yx^T %*% covx^(-1) %*% yx
           = yx^T %*% (covCholx^T %*% covCholx)^(-1) %*% yx
	   = (covCholx^(-T) %*% yx)^T %*% covCholx^(-T) %*% yx
	   = dummy2^T %*% dummy2

     scaleMat = mahal / (k + 1) * dummy3.

  */

  double mahal = 0;
  for (int i=0;i<*ntau;i++)
    mahal += dummy2[i] * dummy2[i];

  mahal /= *DoF;

  double minusMahal = - mahal;

  memcpy(scaleMat, covs, *ntaubar * *ntaubar * sizeof(double));
  F77_CALL(dsyrk)("U", "N", ntaubar, ntau, &minusMahal, dummy1, ntaubar,
		  &mahal, scaleMat, ntaubar
		  FCONE FCONE);// Watch out only the upper part is stored.

  // Fill the lower triangular part
  for (int i=0;i<*ntaubar;i++)
    for (int j=i;j<*ntaubar;j++)
      scaleMat[j + i * *ntaubar] = scaleMat[i + *ntaubar * j];

  free(covCholx); free(covs); free(covsx); free(yx); free(dummy1); free(dummy2);
  return;
}

void condsimschlather(int *nsim, int *n, int *nCond, int *allPart, double *cov,
		      double *y, double *sim, double *subextfct,
		      double *extfct, double *timings){

  clock_t start, end;
  int oneInt = 1,
    *currentPart = malloc(*nCond * sizeof(int));
  double *x = malloc(*n * sizeof(double)),
    *z = malloc(*n * sizeof(double)),
    *prop = malloc(*n * sizeof(double));

  /* Get the Cholesky decomposition of the covariance matrix */
  int info = 0;
  double *covChol = malloc(*n * *n * sizeof(double));
  memcpy(covChol, cov, *n * *n * sizeof(double));
  F77_CALL(dpotrf)("U", n, covChol, n, &info FCONE);

  if (info != 0)
    error("Error code %d from Lapack routine '%s'", info, "dpotrf");

  GetRNGstate();

  for (int i=0;i<*nsim;i++){

    /* Select the partition */
    memcpy(currentPart, allPart + i * *nCond, *nCond * sizeof(int));
    int size = getPartSize(currentPart, nCond);

    start = clock();
    for (int j=0;j<*n;j++)
      z[j] = -1e6;

    for (int j=0;j<size;j++){

      int r = 0;
      for (int k=0;k<*nCond;k++)
	r += (currentPart[k] == j);

      int nr = *n - r,
	*tau = malloc(r * sizeof(int));

      gettau(nCond, currentPart, &j, tau);

      /* Watch out taubar typically contains some conditioning
	 locations and the locations where we want to get conditional
	 realizations!!! */
      int *taubar = malloc(nr * sizeof(int));
      for (int k=0;k<(*n - *nCond);k++)
	taubar[*nCond - r + k] = *nCond + k;

      if (r < *nCond)
	gettaubar(nCond, currentPart, &j, taubar);

      double DoF = 0,
	*scaleMat = malloc(nr * nr * sizeof(double)),
	*mu = malloc(nr * sizeof(double));
      getParametersSC(tau, taubar, &r, &nr, cov, y, &DoF, mu, scaleMat);

      // Compute the Cholesky decomposition of scaleMat
      double *scaleMatChol = malloc(nr * nr * sizeof(double));
      memcpy(scaleMatChol, scaleMat, nr * nr * sizeof(double));

      int info = 0;
      F77_CALL(dpotrf)("U", &nr, scaleMatChol, &nr, &info FCONE);

      if (info != 0)
	error("2. error code %d from Lapack routine '%s'", info, "dpotrf");

      /* Simulate student processes that satisfy the upper bound
	 constraints from iBchol */
      for (int k=0;k<r;k++)
	z[tau[k]] = y[tau[k]];

      double *eps = malloc(nr * sizeof(double));
      int flag = 1;

      while (flag){
	flag = *nCond - r;

	double scaleFactor = sqrt(DoF / rchisq(DoF));

	for (int k=0;k<nr;k++)
	  eps[k] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", &nr, scaleMatChol, &nr, eps, &oneInt
	     FCONE FCONE FCONE);

	for (int k=0;k<nr;k++)
	  eps[k] = mu[k] + scaleFactor * eps[k];

	for (int k=0;k<(*nCond - r);k++)
	  flag -= (eps[k] <= y[taubar[k]]);
      }

      for (int k=0;k<nr;k++)
	z[taubar[k]] = fmax2(z[taubar[k]], eps[k]);

      free(scaleMatChol); free(tau); free(taubar); free(mu); free(scaleMat);
      free(eps);
    }

    end = clock();
    timings[0] = ((double) (end - start)) / CLOCKS_PER_SEC;

    /* Simulate Schlather processes that satisfy the upper bound
       constraints */
    start = clock();
    const double normCst = sqrt(M_2PI);
    double poisson = 0;

    for (int j=0;j<*n;j++)
      x[j] = -1e6;

    int nKO = *n;
    while (nKO){

      poisson += exp_rand();
      double ipoisson = normCst / poisson,
	thresh = 3.5 * M_2PI * ipoisson;

      for (int j=0;j<*n;j++)
	prop[j] = norm_rand();

      F77_CALL(dtrmv)("U", "T", "N", n, covChol, n, prop, &oneInt
		      FCONE FCONE FCONE);

      for (int j=0;j<*n;j++)
	prop[j] *= ipoisson;

      int flag = *nCond;
      for (int j=0;j<*nCond;j++)
	flag -= (prop[j] <= y[j]);

      nKO = *n;
      if (flag == 0){
	for (int j=0;j<*n;j++)
	  x[j] = fmax2(x[j], prop[j]);

	//nKO--;
	for (int j=0;j<*n;j++)
	  nKO -= (thresh <= x[j]);
      }
    }

    end = clock();
    timings[1] = ((double) (end - start)) / CLOCKS_PER_SEC;

    /* Get the max between those two simulations */
    for (int j=0;j<*n;j++){
      sim[i * *n + j] = fmax2(x[j], z[j]);
      subextfct[i * *n + j] = x[j];
      extfct[i * *n + j] = z[j];
    }
  }

  PutRNGstate();

  free(currentPart); free(x); free(z); free(prop); free(covChol);
  return;
}

void gettau(int *nCond, int *part, int *set, int *tau){

  /* This functions get the elements that belong to a given set of a
     partition. */
  int idx1 = 0;
  for (int i=0;i<*nCond;i++){
    if (part[i] == *set){
      tau[idx1] = i;
      idx1++;
    }
  }

  return;
}

void gettaubar(int *nCond, int *part, int *set, int *taubar){

  /* This functions get the elements that do not belong to a given set
     of a partition. */

  int idx1 = 0;
  for (int i=0;i<*nCond;i++){
    if (part[i] != *set){
      taubar[idx1] = i;
      idx1++;
    }
  }

  return;
}

void gibbsForPartSC(int *nchain, int *nthin, int *burnin, int *nCond,
		  int *currentPart, double *cov, double *y, int *chain,
		  double *timings){

  /* This function generates a Markov chain using a Gibbs sampler
     whose target distribution is \Pr[\theata = \tau | Z(x) = z].  The
     update is performed by picking up randomly one conditional
     location, compute all the weights related to moving this location
     to another set (or even a new one), and sample the next state of
     the chain from the conditional distribution. */

  clock_t start, end;
  start = clock();
  int *newPart = malloc(*nCond * sizeof(int)),
    size = getPartSize(currentPart, nCond),
    iterThin = 0,
    iter = 0;

  GetRNGstate();

  while (iterThin < *nchain){

    void R_CheckUserInterrupt(void);

    /* Pick up one location randomly */
    int idxSite = *nCond * unif_rand(), currentSet = currentPart[idxSite];
    memcpy(newPart, currentPart, *nCond * sizeof(int));

    int r = 0;
    for (int j=0;j<*nCond;j++)
      r += (currentPart[j] == currentSet);

    int ncondWeights = size + (r > 1);
    /* Rmk: If r = 1, then we cannot have newSet = size since this
       partition has already been considered when newSet = oldSet */

    double *condWeights = malloc(ncondWeights * sizeof(double));

    for (int newSet=0;newSet<ncondWeights;newSet++){
      if (newSet == currentSet){
	/* This means that the site is moved to the same partition set
	   --> the partition is left unchanged */
	condWeights[newSet] = 1;
	continue;
      }

      newPart[idxSite] = newSet;

      condWeights[newSet] = computeWeightForOneSetSC(&newSet, nCond, newPart, cov, y) /
	computeWeightForOneSetSC(&currentSet, nCond, currentPart, cov, y);

      if (r > 1)
	condWeights[newSet] *= computeWeightForOneSetSC(&currentSet, nCond, newPart, cov, y);

      int r2 = 0;
      for (int k=0;k<*nCond;k++)
	r2 += (currentPart[k] == newSet);

      if (r2 > 0)
	condWeights[newSet] /= computeWeightForOneSetSC(&newSet, nCond, currentPart, cov, y);
    }

    /* Compute the normalizing constant and normalize the
       distribution. */
    double isumCondWeights = 0;
    for (int j=0;j<ncondWeights;j++)
      isumCondWeights += condWeights[j];

    isumCondWeights = 1 / isumCondWeights;

    for (int j=0;j<ncondWeights;j++)
      condWeights[j] *= isumCondWeights;

    /* Sample from the full conditional distribution. */
    double u = unif_rand();
    int newSet = -1;
    {
      int flag = 1;
      double cumWeight = 0;
      while (flag){
	newSet++;
	cumWeight += condWeights[newSet];
	flag = (u >= cumWeight);
      }
    }

    /* Update the Markov chain. */
    if (newSet != currentSet){
      currentPart[idxSite] = newSet;

      if (r == 1)
	size--;

      else
	size += (newSet == size);

      convert2rightformat(currentPart, nCond, &size);
    }

    if ((iter > *burnin) & ((iter % *nthin) == 0)){
      memcpy(chain + iterThin * *nCond, currentPart, *nCond * sizeof(int));
      iterThin++;
    }

    iter++;
    free(condWeights);
  }

  PutRNGstate();

  free(newPart);

  end = clock();
  *timings = ((double) (end - start)) / CLOCKS_PER_SEC;

  return;
}

double computeWeightForOneSetSC(int *idxSet, int *nCond, int *partition,
			      double *cov, double *y){

  /* This function computes the weights of a single set of a given
     partition. Hence the weight of a partition is the product of
     these single weights. (Schlather case) */
  double weight = 1;
  int r = 0;

  for (int i=0;i<*nCond;i++)
    r += (partition[i] == *idxSet);

  int nr = *nCond - r,
    *tau = malloc(r * sizeof(int)),
    *taubar = malloc(nr * sizeof(int));

  gettau(nCond, partition, idxSet, tau);

  if (r < *nCond) {
    gettaubar(nCond, partition, idxSet, taubar);

    /* Define the parameters of the Student distribution, i.e., the
       degree of freedom (DoF), the Cholesky decomposition of the
       (scaleMat) scale matrix and the location parameter (mu) */
    double DoF = 0,
      *scaleMat = malloc(nr * nr * sizeof(double)),
      *mu = malloc(nr * sizeof(double));
    getParametersSC(tau, taubar, &r, &nr, cov, y, &DoF, mu, scaleMat);

    /* Compute the proba values --- see the note */
    double prob = 0;
    computeprobaSC(&DoF, mu, scaleMat, y, &nr, taubar, &prob);
    weight = prob;

    free(scaleMat); free(mu);
  }

  /* Compute the f values --- see the note */
  double f = 0;
  getfvaluesSC(y, nCond, &r, tau, cov, &f);
  // f values are OK

  weight *= exp(f);

  free(tau); free(taubar);
  return weight;
}

void convert2rightformat(int *partition, int *n, int *size){

  int bound = 0;

  for (int i=0;i<*n-1;i++){

    if (partition[i] > bound){
      int old = partition[i], new = bound;

      for (int j=i;j<*n;j++){

	if (partition[j] == old){
	  partition[j] = new;
	  continue;
	}

	if (partition[j] == new)
	  partition[j] = old;

      }
    }

    bound = 0;
    for (int j=0;j<=i;j++)
      bound = imax2(bound, partition[j]);

    bound++;
  }

  /* The following lines work if we know the partition size

     if (partition[*n-1] >= *size)
     partition[*n-1] = *size - 1;

  */

  // Temptative to avoid knowing the partition size
  if (partition[*n-1] > bound)
    partition[*n-1] = bound;

  return;
}


void sampleDiscreteDist(int *n, double *prob, int *ans){

  GetRNGstate();
  double u = unif_rand();

  int idx = -1;
  {
    int flag = 1;
    double bound = 0;
    while (flag){
      idx++;
      bound += prob[idx];
      flag = (u >= bound);
    }
  }

  PutRNGstate();

  *ans = idx;
  return;
}

void validPart(int *partition, int *n, int *valid){

  int bound = 0;
  *valid = 1;

  for (int i=0;i<*n-1;i++){

    if (partition[i] > bound){
      *valid = 0;
      return;
    }

    bound = 0;
    for (int j=0;j<=i;j++)
      bound = imax2(bound, partition[j]);

    bound++;
  }

  if (partition[*n-1] > bound)
    *valid = 0;

  return;
}

void getBounds(int *partition, int *n, int *idx, int *lbound, int *ubound){
  /* This function computes the lower and upper bounds for the values
     that are admissible for partition[idx] to be a valid partition */

  *ubound = *lbound = 0;
  for (int i=1;i<*idx;i++)
    *ubound = imax2(*ubound, partition[i]);

  *ubound = *ubound + 1;

  if (*idx < (*n-1)) {
    for (int i=*idx+1;i<*n;i++){
      if (partition[i] == *ubound)
	break;

      else if (partition[i] > *ubound){
	*lbound = partition[i] - 1;
	break;
      }
    }
  }

  return;
}

int getPartSize(int *partition, int *n){

  /* This function computes the size of a partition */
  int size = 0;

  for (int i=0;i<*n;i++)
    size = imax2(size, partition[i]);

  size++;

  return size;
}

void totoSC(int *nsim, int *n, double *covChol, double *ans){

  int oneInt = 1;
  const double normCst = M_SQRT2 * M_SQRT_PI;

  double *x = malloc(*n * sizeof(double)),
    *prop = malloc(*n * sizeof(double));

  for (int i=0;i<*nsim;i++){
    for (int j=0;j<*n;j++)
      x[j] = 0;

    double poisson = 0;
    int nKO = *n;

    while(nKO){
      poisson += exp_rand();

      double ipoisson = 1 / poisson,
	thresh = 3.5 * normCst * ipoisson;

      for (int j=0;j<*n;j++)
	prop[j] = norm_rand();

      F77_CALL(dtrmv)("U", "T", "N", n, covChol, n, prop, &oneInt FCONE FCONE FCONE);

      for (int j=0;j<*n;j++)
	prop[j] = fmax2(0, normCst * prop[j] * ipoisson);

      for (int j=0;j<*n;j++)
	x[j] = fmax2(x[j], prop[j]);

      nKO = *n;
      for (int j=0;j<*n;j++)
	nKO -= (thresh <= x[j]);
    }

    for (int j=0;j<*n;j++)
      ans[i + j * *nsim] = x[j];
  }

  free(x); free(prop);

  return;
}

void getStartingPartitionSC(int *nsim, int *n, double *covChol, int *startPart){

  /* This function simulates a bunch of max-stable realizations and
     define a suitable starting values for the Gibbs sampler. */
  int oneInt = 1;
  const double normCst = M_SQRT2 * M_SQRT_PI;

  int *dummy = malloc(*n * sizeof(int));
  double *x = malloc(*n * sizeof(double)),
    *prop = malloc(*n * sizeof(double));

  for (int i=0;i<*nsim;i++){
    for (int j=0;j<*n;j++){
      x[j] = 0;
      dummy[j] = -1;
    }

    double poisson = 0;
    int nKO = *n;

    int partSet = 0;
    while(nKO){
      poisson += exp_rand();

      double ipoisson = 1 / poisson,
	thresh = 3.5 * normCst * ipoisson;

      for (int j=0;j<*n;j++)
	prop[j] = norm_rand();

      F77_CALL(dtrmv)("U", "T", "N", n, covChol, n, prop, &oneInt FCONE FCONE FCONE);

      for (int j=0;j<*n;j++)
	prop[j] = fmax2(0, normCst * prop[j] * ipoisson);

      int hasChanged = 0;
      for (int j=0;j<*n;j++){
	if (prop[j] > x[j]){
	  dummy[j] = partSet;
	  hasChanged = 1;
	}

	x[j] = fmax2(x[j], prop[j]);
      }

      nKO = *n;
      for (int j=0;j<*n;j++)
	nKO -= (thresh <= x[j]);

      if (hasChanged){
	partSet++;
	convert2rightformat(dummy, n, &partSet);
	partSet = getPartSize(dummy, n);
      }
    }

    for (int j=0;j<*n;j++)
      startPart[i * *n + j] = dummy[j];
  }

  free(x); free(prop); free(dummy);

  return;
}

void pmvnorm(double *bounds, int *n, double *cor, double *prob,
	     double *err, int *nMc){
  /* This function computes the CDF of a standard multivariate normal
     random vector using quasi monte carlo sequence*/

  const int primeNumbers[100] = {
    2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
    31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
    73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
    353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
    419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
    467, 479, 487, 491, 499, 503, 509, 521, 523, 541};

  // 1. We need to reorder the variables
  double *newCor = malloc(*n * *n * sizeof(double));
  memcpy(newCor, cor, *n * *n * sizeof(double));

  int *idx = malloc(*n * sizeof(int));
  for (int i=0;i<*n;i++)
    idx[i] = i;

  rsort_with_index(bounds, idx, *n);

  for (int i=0;i<*n;i++)
    for (int j=0;j<*n;j++)
      newCor[i + j * *n] = cor[idx[i] + *n * idx[j]];

  //Compute the Cholesky decomposition of newCor
  int info = 0;
  F77_CALL(dpotrf)("L", n, newCor, n, &info FCONE);

  if (info != 0)
    error("1. error code %d from Lapack routine '%s'", info, "dpotrf");

  // Get the inverse entries of newCor for efficiency
  double *inewCor = malloc(*n * *n * sizeof(double));
  for (int i=0;i<*n * *n;i++)
    inewCor[i] = 1 / newCor[i];

  // Initialize objects
  int M = 1, N = *n * 1000;
  double T = 0, E = 1, V = 0;
  double *I = malloc(M * sizeof(double)),
    *Delta = malloc(*n * sizeof(double)),
    delta,
    *q = malloc(*n * sizeof(double)),
    *e = malloc(*n * sizeof(double)),
    *f = malloc(*n * sizeof(double)),
    *w = malloc(*n * sizeof(double)),
    *y = malloc((*n - 1) * sizeof(double)),
    *ea = malloc(*n * sizeof(double)),
    *fa = malloc(*n * sizeof(double)),
    *wa = malloc(*n * sizeof(double)),
    *ya = malloc((*n - 1) * sizeof(double));//a for antithethic

  for (int i=0;i<*n;i++)
    q[i] = sqrt(primeNumbers[i]);

  e[0] = ea[0] = pnorm(bounds[0] * inewCor[0], 0, 1, 1, 0);
  f[0] = fa[0] = e[0];

  GetRNGstate();

  for (int i=0;i<M;i++){
    I[i] = 0;

    for (int j=0;j<*n;j++)
      Delta[j] = unif_rand();

    for (int j=0;j<N;j++){
      for (int k=0;k<*n;k++){
	double dummy = j * q[k] + Delta[k];
	w[k] = fabs(2 * (dummy - ftrunc(dummy)) - 1);
	wa[k] = 1 - w[k];
      }

      for (int k=1;k<*n;k++){
	y[k-1] = qnorm(w[k-1] * e[k-1], 0, 1, 1, 0);
	ya[k-1] = qnorm(wa[k-1] * ea[k-1], 0, 1, 1, 0);

	double dummy = 0, dummya = 0;
	for (int l=0;l<k;l++){
	  dummy += newCor[k + l * *n] * y[l];
	  dummya += newCor[k + l * *n] * ya[l];
	}

	e[k] = pnorm((bounds[k] - dummy) * inewCor[k * (*n + 1)], 0, 1, 1, 0);
	ea[k] = pnorm((bounds[k] - dummya) * inewCor[k * (*n + 1)], 0, 1, 1, 0);
	f[k] = e[k] * f[k - 1];
	fa[k] = ea[k] * fa[k - 1];
      }

      I[i] = I[i] + (0.5 * (f[*n - 1] + fa[*n - 1]) - I[i]) / ((double) j + 1);
    }

    delta = (I[i] - T) / ((double) i + 1);
    T += delta;
    V = (i - 1) * V / ((double) i + 1) + delta * delta;
    E = 3 * sqrt(V);
  }

  PutRNGstate();

  *prob = T;
  *err = E;
  *nMc = (int )N;

  free(newCor); free(idx); free(inewCor); free(I); free(Delta);
  free(q); free(e); free(f); free(w); free(y); free(ea); free(fa);
  free(wa); free(ya);

  return;
}

void computeWeightsBR(int *nCond, double *y, int *nPart, int *allPart, int *allSize,
		      double *cov, double *sigma2, double *covChol, double *ham,
		      double *mean1, double *weights){
  /* This function computes the distribution function of the random
     partition given the conditioning values, i.e., Pr[\theta = \tau |
     Z(x) = z]. (Brown-Resnick case) */

  int info = 0;
  int *currentPart = malloc(*nCond *  sizeof(int));

  // Set the weights to one
  for (int i=0;i<*nPart;i++)
    weights[i] = 1;

  for (int i=0;i<*nPart;i++){
    // Loop over all the partitions
    memcpy(currentPart, allPart + i * *nCond, *nCond * sizeof(int));

    for (int j=0;j<allSize[i];j++){
      //Loop over all sets of the current partition

      int r = 0;
      for (int k=0;k<*nCond;k++)
	r += (currentPart[k] == j);

      int nr = *nCond - r,
	*tau = malloc(r * sizeof(int)),
	*taubar = malloc(nr * sizeof(int));

      gettau(nCond, currentPart, &j, tau);

      if (r < *nCond) {
	gettaubar(nCond, currentPart, &j, taubar);

	// Build the matrices J and Jtilde
	double *J = malloc(*nCond * nr * sizeof(double)),
	  *Jtilde = malloc(*nCond * *nCond * sizeof(double));
	buildJ(tau, nCond, &r, J);
	buildJtilde(tau, nCond, &r, Jtilde);

	/* Define the parameters of the distribution of the extremal
	   function. For the Brown-Resnick case, it is a multivariate
	   log-normal distribution. Hence below we compute the
	   Cholesky decomposition of the **INVERSE** of the covariance
	   matrix and the mean mu */
	double *iBchol = malloc(nr * nr * sizeof(double)),
	  *mu = malloc(nr * sizeof(double));
	getParametersBR(J, Jtilde, nCond, &nr, covChol, ham, mean1, y, iBchol,
			mu);

	/* Compute the proba values, i.e., the probabilities that one
	   extremal function is smaller than the hitting bounds at the
	   remaining conditioning locations. */
	double prob = 0;
	computeprobaBR(iBchol, mu, y, nCond, &r, nCond, taubar, &prob);
	weights[i] *= prob;

	free(J); free(Jtilde); free(iBchol); free(mu);
      }

      /* Get the sub-matrix of cov */
      double *covjchol = malloc(r * r * sizeof(double));
      getSubMatrix(cov, nCond, &r, tau, &r, tau, covjchol);

      /* Get the Cholesky decomposition of this submatrix */
      F77_CALL(dpotrf)("U", &r, covjchol, &r, &info FCONE);

      if (info != 0)
	error("4. error code %d from Lapack routine '%s'", info, "dpotrf");

      //Get the sub-vector of sigma2 and of y
      double *sigma2j = malloc(r * sizeof(double)),
	*yj = malloc(r * sizeof(double));
      for (int k=0;k<r;k++){
	sigma2j[k] = sigma2[tau[k]];
	yj[k] = y[tau[k]];
      }

      /* Compute the intensity function values, i.e.,
	 \lambda_{x_{\tau_j}}(z_{\tau_j}). */
      double f = 0;
      getfvaluesBR(yj, sigma2j, covjchol, &r, &f);

      weights[i] *= exp(f);
      free(tau); free(taubar); free(sigma2j); free(yj);
    }
  }

  /* Compute the normalizing constant and normalize the
     distribution */
  double isumWeights = 0;
  for (int i=0;i<*nPart;i++)
    isumWeights += weights[i];

  isumWeights = 1 / isumWeights;

  for (int i=0;i<*nPart;i++)
    weights[i] *= isumWeights;

  free(currentPart);
  return;
}

void computeprobaBR(double *icovChol, double *mu, double *y, int *n,
		    int *r, int *nCond, int *taubar, double *prob){

  int nr = *n - *r, ncondr = *nCond - *r;
  *prob = 0;

  /* So far we only have compute the inverse of B. We need to inverse
     it right now using iBchol */

  double *B = malloc(nr * nr * sizeof(double));
  memcpy(B, icovChol, nr * nr * sizeof(double));

  {
    int info = 0;
    F77_CALL(dpotri)("U", &nr, B, &nr, &info FCONE);

    if (info != 0)
      error("5. error code %d from Lapack routine '%s'", info, "dpotri");
  }

  /* Get the submatrix of B */
  double *Bj = malloc(ncondr * ncondr * sizeof(double));
  int *idx = malloc(ncondr * sizeof(int));
  for (int i=0;i<ncondr;i++)
    idx[i] = i;

  getSubMatrix(B, &nr, &ncondr, idx, &ncondr, idx, Bj);

  /* Since Bj is symmetric we need to fill in the lower diagonal
     elements */
  for (int i=0;i<ncondr;i++)
    for (int j=0;j<i;j++)
      Bj[i + ncondr * j] = Bj[j + ncondr * i];

  /* Get the subvector of mu and y*/
  double *muj = malloc(ncondr * sizeof(double)),
    *yj = malloc(ncondr * sizeof(double));

  for (int i=0;i<ncondr;i++){
    muj[i] = mu[i];
    yj[i] = y[taubar[i]];
  }

  if (ncondr == 1)
    *prob = pnorm(yj[0], muj[0], sqrt(Bj[0]), 1, 0);

  else {
    /* Standardize Bj and yj before the call to pmvnorm */
    standardize(yj, Bj, muj, &ncondr);

    /* Computing the proba */
    double err = 0;
    int nMc = 0;

    pmvnorm(yj, &ncondr, Bj, prob, &err, &nMc);
  }

  free(B); free(Bj); free(idx); free(muj); free(yj);

  return;
}

void standardize(double *quant, double *cov, double *mean, int *n){
  /* Standardize the covariance matrix first */
  double *isd = malloc(*n * sizeof(double));

  for (int i=0;i<*n;i++)
    isd[i] = 1 / sqrt(cov[i * (1 + *n)]);

  for (int i=0;i<*n;i++)
    for (int j=0;j<*n;j++)
      cov[i + *n * j] *= isd[i] * isd[j];

  /* Standardize quant */
  for (int i=0;i<*n;i++)
    quant[i] = (quant[i] - mean[i]) * isd[i];

  free(isd);
}

void getfvaluesBR(double *y, double *sigma2, double *covjchol, int *r,
		  double *f){

  /* Compute the log of the determinant of covj */
  int oneInt = 1;
  double detcovj = 0;
  for (int i=0;i<*r;i++)
    detcovj += log(covjchol[i * (1 + *r)]);

  detcovj *= 2;

  double *mean = malloc(*r * sizeof(double));
  for (int i=0;i<*r;i++)
    mean[i] = y[i] + 0.5 * sigma2[i];

  double *one = malloc(*r * sizeof(double));
  for (int i=0;i<*r;i++)
    one[i] = 1;

  /* Compute covjchol^{-T} %*% 1 and covjchol^{-T} %*% mean*/
  double *covjchol1 = malloc(*r * sizeof(double)),
    *covjcholmean = malloc(*r * sizeof(double));

  memcpy(covjchol1, one, *r * sizeof(double));
  memcpy(covjcholmean, mean, *r * sizeof(double));

  F77_CALL(dtrsv)("U", "T", "N", r, covjchol, r, covjchol1, &oneInt FCONE FCONE FCONE);
  F77_CALL(dtrsv)("U", "T", "N", r, covjchol, r, covjcholmean, &oneInt FCONE FCONE FCONE);

  double mahal1 = 0, mahalmean = 0, mahal1mean = 0;
  for (int i=0;i<*r;i++){
    mahal1 += covjchol1[i] * covjchol1[i];
    mahalmean += covjcholmean[i] * covjcholmean[i];
    mahal1mean += covjchol1[i] * covjcholmean[i];
  }

  *f = (1 - *r) * M_LN_SQRT_2PI - 0.5 * (detcovj + log(mahal1) + mahalmean -
					 (mahal1mean - 1) * (mahal1mean - 1) / mahal1);

  for (int i=0;i<*r;i++)
    *f -= y[i];

  free(mean); free(one); free(covjchol1); free(covjcholmean);

  return;
}

void buildJ(int *tau, int *n, int *r, double *J){
  /* This function computes the J matrix of Dombry et al. (2012)  */
  for (int k=0;k<(*n * (*n - *r));k++)
    J[k] = 0;

  // Could be optimized I think
  int idx = 0;
  for (int i=0;i<*n;i++){
    int flag = 1;

    for (int j=0;j<*r;j++)
      flag &= (i != tau[j]);

    if (flag){
      J[i + idx * *n] = 1;
      idx++;
    }
  }

  return;
}

void buildJtilde(int *tau, int *n, int *r, double *Jtilde){
  /* This function computes the Jtilde matrix of Dombry et
     al. (2012) */
  for (int i=0;i<(*n * *n);i++)
    Jtilde[i] = 0;

  for (int i=0;i<*r;i++)
    Jtilde[tau[i] * (*n + 1)] = 1;

  return;
}

void getParametersBR(double *J, double *Jtilde, int *n, int *nr,
		     double *covChol, double *ham, double *mean1,
		     double *ytilde, double *iBchol, double *mu){

  /* This function computes the parameters of distribution of the
     extremal functions. For the Brown--Resnick case, it is a multivariate
     log-normal distribution with mean (mu) and inverse covariance matrix (iB). */
  int oneInt = 1, info = 0;
  double alpha = 1, beta = 0;
  /* Compute bread = icovChol^T %*% J, i.e. solve covChol^T %*% X
     = J */
  double *bread = malloc(*n * *nr * sizeof(double));
  memcpy(bread, J, *n * *nr * sizeof(double));

  alpha = 1;
  F77_CALL(dtrsm)("L", "U", "T", "N", n, nr, &alpha, covChol, n,
		  bread, n 
		  FCONE FCONE FCONE FCONE);

  /* Finally the inverse of the covariance matrix is iB = bread^T
     %*% ham %*% bread */
  double *iB = malloc(*nr * *nr * sizeof(double));
  for (int k=0;k<*nr * *nr;k++)
    iB[k] = 0;

  alpha = 1;
  beta = 0;
  {
    double *dummy = malloc(*n * *nr * sizeof(double));
    for (int k=0;k<(*n * *nr);k++)
      dummy[k] = 0;

    // Since ham is symmetric we first compute ham %*% bread
    F77_CALL(dsymm)("L", "U", n, nr, &alpha, ham, n, bread, n, &beta,
		    dummy, n FCONE FCONE);
    // And then bread^T %*% dummy
    F77_CALL(dgemm)("T", "N", nr, nr, n, &alpha, bread, n, dummy,
		    n, &beta, iB, nr  FCONE FCONE);

    free(dummy);
  }

  // Keep only the upper part of the covariance matrix
  for (int k=0;k<*nr;k++)
    for (int l=0;l<k;l++)
      iB[k + *nr * l] = 0;

  /* Compute breadr = icovChol^T %*% Jtilde, i.e. solve covChol^T
     X = Jtilde */
  double *breadr = malloc(*n * *n * sizeof(double));
  memcpy(breadr, Jtilde, *n * *n * sizeof(double));
  alpha = 1;
  F77_CALL(dtrsm)("L", "U", "T", "N", n, n, &alpha, covChol, n,
		  breadr, n FCONE FCONE FCONE FCONE);

  /* Compute the Cholesky decomposition of iB */
  memcpy(iBchol, iB, *nr * *nr * sizeof(double));
  F77_CALL(dpotrf)("U", nr, iBchol, nr, &info FCONE);

  if (info != 0)
    error("3. error code %d from Lapack routine '%s'", info, "dpotrf");

  /* Compute breadl = bread %*% B we do this from the Cholesky
     decomposition of iB */
  double *breadl = malloc(*n * *nr * sizeof(double));
  memcpy(breadl, bread, *n * *nr * sizeof(double));

  // First solve Y iBchol = bread, i.e., Y = bread %*% iBchol^{-1}
  alpha = 1;
  F77_CALL(dtrsm)("R", "U", "N", "N", n, nr, &alpha, iBchol, nr,
		  breadl, n
		  FCONE FCONE FCONE FCONE);
  /* Then solve X iBchol^T = Y, i.e., X = Y %*% iBchol^{-T} =
     bread %*% iBchol^{-1} %*% iBchol^{-T} = bread %*% (iBchol^T %*%
     iBchol)^{-1} = bread %*% iB^{-1} = bread %*% B */

  F77_CALL(dtrsm)("R", "U", "T", "N", n, nr, &alpha, iBchol, nr,
		  breadl, n
		  FCONE FCONE FCONE FCONE);

  /* Finally mu = breadl^T %*% (mean1 - ham %*% breadr %*% ytilde) */
  for (int k=0;k<*nr;k++)
    mu[k] = 0;

  {
    /* We start by computing dummy1 = breadr %*% ytilde */
    double *dummy1 = malloc(*n * sizeof(double));
    for (int k=0;k<*n;k++)
      dummy1[k] = 0;

    alpha = 1;
    beta = 0;
    F77_CALL(dgemv)("N", n, n, &alpha, breadr, n, ytilde, &oneInt, &beta,
		    dummy1, &oneInt FCONE);

    /* Then we compute dummy2 = mean1 - ham %*% dummy1 */
    double *dummy2 = malloc(*n * sizeof(double));
    memcpy(dummy2, mean1, *n * sizeof(double));

    alpha = -1;
    beta = 1;
    F77_CALL(dsymv)("U", n, &alpha, ham, n, dummy1, &oneInt, &beta,
		    dummy2, &oneInt FCONE);

    /* Lastly we compute mu = breadl^T %*% dummy2 */
    alpha = 1;
    beta = 0;
    F77_CALL(dgemv)("T", n, nr, &alpha, breadl, n, dummy2, &oneInt,
		    &beta, mu, &oneInt FCONE);

    free(dummy1); free(dummy2);
  }

  free(iB); free(bread); free(breadr); free(breadl);

  return;
}

void condsimbrown(int *nsim, int *n, int *nCond, int *allPart, double *covChol,
		  double *sigma2, double *ham, double *mean1, double *ytilde,
		  double *coord, double *range, double *smooth, int *dim,
		  double *sim, double *subextfct, double *extfct, double *timings){

  clock_t start, end;
  int oneInt = 1,
    *currentPart = malloc(*nCond * sizeof(int)),
    blockSize = 500;
  double *x = malloc(*n * sizeof(double)),
    *z = malloc(*n * sizeof(double)),
    *prop = malloc(*n * sizeof(double)),
    *vario = malloc(*n * sizeof(double));

  for (int i=0;i<*n;i++){
    vario[i] = 0;
    for (int j=0;j<*dim;j++)
      vario[i] += coord[i + j * *n] * coord[i + j * *n];

    vario[i] = R_pow(sqrt(vario[i]) / *range, *smooth);
  }

  GetRNGstate();

  for (int i=0;i<*nsim;i++){

    /* Select the partition */
    memcpy(currentPart, allPart + i * *nCond, *nCond * sizeof(int));
    int size = getPartSize(currentPart, nCond);

    start = clock();
    for (int j=0;j<*n;j++)
      z[j] = -1e6;

    for (int j=0;j<size;j++){

      int r = 0;
      for (int k=0;k<*nCond;k++)
	r += (currentPart[k] == j);

      int nr = *n - r,
	*tau = malloc(r * sizeof(int)),
	*taubar = malloc((*nCond - r) * sizeof(int));
      gettau(nCond, currentPart, &j, tau);

      if (r < *nCond)
	gettaubar(nCond, currentPart, &j, taubar);

      double *J = malloc(*n * nr * sizeof(double)),
	*Jtilde = malloc(*n * *n * sizeof(double));

      buildJ(tau, n, &r, J);
      buildJtilde(tau, n, &r, Jtilde);

      double *iBchol = malloc(nr * nr * sizeof(double)),
	*mu = malloc(nr * sizeof(double));
      getParametersBR(J, Jtilde, n, &nr, covChol, ham, mean1, ytilde, iBchol,
		      mu);

      /* Simulate log-normal processes that satisfy the upper bound
	 constraints from iBchol !!! */
      for (int k=0;k<r;k++)
	z[tau[k]] = ytilde[tau[k]];

      double *eps = malloc(nr * sizeof(double));
      int flag = 1;

      while (flag){
	flag = *nCond - r;
	for (int k=0;k<nr;k++)
	  eps[k] = norm_rand();

	F77_CALL(dtrsv)("U", "N", "N", &nr, iBchol, &nr, eps, &oneInt
	     FCONE FCONE FCONE);

	for (int k=0;k<nr;k++)
	  eps[k] += mu[k];

	for (int k=0;k<(*nCond - r);k++)
	  flag -= (eps[k] <= ytilde[taubar[k]]);
      }

      for (int k=0;k<(*nCond - r);k++)
	z[taubar[k]] = fmax2(z[taubar[k]], eps[k]);

      for (int k=*nCond;k<*n;k++)
	z[k] = fmax2(z[k], eps[k - r]);

      free(tau); free(taubar); free(J); free(Jtilde); free(iBchol); free(mu); free(eps);
    }

    end = clock();
    timings[0] = ((double) (end - start)) / CLOCKS_PER_SEC;

    /* Simulate Brown--Resnick processes that satisfy the upper bound
       constraints */
    start = clock();
    for (int j=0;j<*n;j++)
      x[j] = -1e6;


    double poisson = 0;
    int nBlock = blockSize;

    while(nBlock){
      poisson += exp_rand();

      for (int j=0;j<*n;j++)
	prop[j] = norm_rand();

      F77_CALL(dtrmv)("U", "T", "N", n, covChol, n, prop, &oneInt
		      FCONE FCONE FCONE);

      for (int j=0;j<*n;j++)
	prop[j] = prop[j] - vario[j] - log(poisson);

      int flag = *nCond;
      for (int j=0;j<*nCond;j++)
	flag -= (prop[j] <= ytilde[j]);

      if (flag == 0){
	for (int j=0;j<*n;j++)
	  x[j] = fmax2(x[j], prop[j]);

	nBlock--;
      }
    }

    end = clock();
    timings[1] = ((double) (end - start)) / CLOCKS_PER_SEC;

    /* Get the max between those two simulations and go back to the
       Frechet scale*/
    for (int j=0;j<*n;j++){
      sim[i * *n + j] = exp(fmax2(x[j], z[j]));
      subextfct[i * *n + j] = exp(x[j]);
      extfct[i * *n + j] = exp(z[j]);
    }

  }

  PutRNGstate();

  free(currentPart); free(x); free(z); free(prop); free(vario);

  return;
}

void condsimbrown2(int *nsim, int *n, int *nCond, int *allPart, double *covChol,
		   double *sigma2, double *ham, double *mean1, double *ytilde,
		   double *sim, double *coord, double *range, double *smooth,
		   double *xlim){

  int oneInt = 1, currentPart[*nCond], blockSize = 100, nPoissSup = 200;
  double x[*n], z[*n], prop[*n], varioShift[*n], bounds[2],
    irange = 1 / *range, inflate = *range * R_pow(2, 1 / *smooth),
    covShift[*n * *n];

  bounds[0] = xlim[0] - inflate;
  bounds[1] = xlim[1] + inflate;

  GetRNGstate();
  for (int i=0;i<*nsim;i++){

    /* Select the partition */
    memcpy(currentPart, allPart + i * *nCond, *nCond * sizeof(int));
    int size = getPartSize(currentPart, nCond);

    for (int j=0;j<*n;j++)
      z[j] = -1e6;

    for (int j=0;j<size;j++){

      int r = 0;
      for (int k=0;k<*nCond;k++)
	r += (currentPart[k] == j);

      int nr = *n - r, tau[r], taubar[*nCond - r];
      gettau(nCond, currentPart, &j, tau);

      if (r < *nCond)
	gettaubar(nCond, currentPart, &j, taubar);

      double J[*n * nr], Jtilde[*n * *n];

      buildJ(tau, n, &r, J);
      buildJtilde(tau, n, &r, Jtilde);

      double iBchol[nr * nr], mu[nr];
      getParametersBR(J, Jtilde, n, &nr, covChol, ham, mean1, ytilde, iBchol,
		      mu);

      /* Simulate log-normal processes that satisfy the upper bound
	 constraints from iBchol !!! */
      for (int k=0;k<r;k++)
	z[tau[k]] = ytilde[tau[k]];

      double eps[nr];
      int flag = 1;

      while (flag){
	flag = *nCond - r;
	for (int k=0;k<nr;k++)
	  eps[k] = norm_rand();

	F77_CALL(dtrsv)("U", "N", "N", &nr, iBchol, &nr, eps, &oneInt
			FCONE FCONE FCONE);

	for (int k=0;k<nr;k++)
	  eps[k] += mu[k];


	for (int k=0;k<(*nCond - r);k++)
	  flag -= (eps[k] <= ytilde[taubar[k]]);
      }

      for (int k=0;k<(*nCond - r);k++)
	z[taubar[k]] = fmax2(z[taubar[k]], eps[k]);

      for (int k=*nCond;k<*n;k++)
	z[k] = fmax2(z[k], eps[k - r]);

    }

    /* Simulate Brown--Resnick processes that satisfy the upper bound
       constraints */

    for (int j=0;j<*n;j++)
      x[j] = -1e6;

    for (int nSup=0;nSup<nPoissSup;nSup++){
      double poisson = 0;
      int nBlock = blockSize;

      while(nBlock){

	double xShift = bounds[0] + (bounds[1]- bounds[0]) * unif_rand();

	for (int j=0;j<*n;j++)
	  varioShift[j] = R_pow(fabs(coord[j] - xShift) * irange, *smooth);

	for (int j=0;j<*n;j++){
	  for (int k=0;k<=j;k++){
	    covShift[j + *n * k] = 0;
	    covShift[k + *n * j] = varioShift[j] + varioShift[k] -
	      R_pow(fabs(coord[j] - coord[k]) * irange, *smooth);
	  }
	}

	//Cholesky decomposition
	int info = 0;
	F77_CALL(dpotrf)("U", n, covShift, n, &info FCONE);

	if (info != 0)
	  error("Error code %d from Lapack routine '%s'", info, "dpotrf");

	for (int j=0;j<*n;j++)
	  prop[j] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", n, covShift, n, prop, &oneInt
	     FCONE FCONE FCONE);

	poisson += nPoissSup * exp_rand();
	for (int j=0;j<*n;j++)
	  prop[j] = prop[j] - varioShift[j] - log(poisson);

	int flag = *nCond;
	for (int j=0;j<*nCond;j++)
	  flag -= (prop[j] <= ytilde[j]);

	if (flag == 0){
	  for (int j=0;j<*n;j++)
	    x[j] = fmax2(x[j], prop[j]);

	  nBlock--;
	}
      }
    }

    /* Get the max between those two simulations and go back to the
       Frechet scale*/
    for (int j=0;j<*n;j++)
      sim[i * *n + j] = exp(fmax2(x[j], z[j]));

  }

  PutRNGstate();

  return;
}

void gibbsForPartBR(int *nchain, int *nthin, int *burnin, int *nCond,
		    int *currentPart, double *cov, double *sigma2,
		    double *covChol, double *ham, double *mean1, double *y,
		    int *chain, double *timing){

  /* This function generates a Markov chain using a Gibbs sampler
     whose target distribution is \Pr[\theata = \tau | Z(x) = z].  The
     update is performed by picking up randomly one conditional
     location, compute all the weights related to moving this location
     to another set (or even a new one), and sample the next state of
     the chain from the conditional distribution. */

  clock_t start, end;
  start = clock();

  int *newPart = malloc(*nCond * sizeof(int)),
    size = getPartSize(currentPart, nCond),
    iterThin = 0,
    iter = 0;

  GetRNGstate();

  while (iterThin < *nchain){

    void R_CheckUserInterrupt(void);

    /* Pick up one location randomly */
    int idxSite = *nCond * unif_rand(), currentSet = currentPart[idxSite];
    memcpy(newPart, currentPart, *nCond * sizeof(int));

    int r = 0;
    for (int j=0;j<*nCond;j++)
      r += (currentPart[j] == currentSet);

    int ncondWeights = size + (r > 1);
    /* Rmk: If r = 1, then we cannot have newSet = size since this
       partition has already been considered when newSet = oldSet */

    double *condWeights = malloc(ncondWeights * sizeof(double));

    for (int newSet=0;newSet<ncondWeights;newSet++){

      if (newSet == currentSet){
	/* This means that the site is moved to the same partition set
	   --> the partition is left unchanged */
	condWeights[newSet] = 1;
	continue;
      }

      newPart[idxSite] = newSet;

      condWeights[newSet] = computeWeightForOneSetBR(&newSet, nCond, newPart, cov, sigma2,
						     covChol, ham, mean1, y) /
	computeWeightForOneSetBR(&currentSet, nCond, currentPart, cov, sigma2,
				 covChol, ham, mean1, y);

      if (r > 1)
	condWeights[newSet] *= computeWeightForOneSetBR(&currentSet, nCond, newPart, cov, sigma2,
							covChol, ham, mean1, y);

      int r2 = 0;
      for (int k=0;k<*nCond;k++)
	r2 += (currentPart[k] == newSet);

      if (r2 > 0)
	condWeights[newSet] /= computeWeightForOneSetBR(&newSet, nCond, currentPart, cov, sigma2,
							covChol, ham, mean1, y);
    }

    /* Compute the normalizing constant and normalize the
       distribution. */
    double isumCondWeights = 0;
    for (int j=0;j<ncondWeights;j++)
      isumCondWeights += condWeights[j];

    isumCondWeights = 1 / isumCondWeights;

    for (int j=0;j<ncondWeights;j++)
      condWeights[j] *= isumCondWeights;

    /* Sample from the full conditional distribution. */
    double u = unif_rand();

    int newSet = -1;
    {
      int flag = 1;
      double cumWeight = 0;
      while (flag){
	newSet++;
	cumWeight += condWeights[newSet];
	flag = (u >= cumWeight);
      }
    }

    /* Update the Markov chain. */

    if (newSet != currentSet){
      currentPart[idxSite] = newSet;

      if (r == 1)
	size--;

      else
	size += (newSet == size);

      convert2rightformat(currentPart, nCond, &size);
    }

    if ((iter > *burnin) & ((iter % *nthin) == 0)){
      memcpy(chain + iterThin * *nCond, currentPart, *nCond * sizeof(int));
      iterThin++;
    }

    iter++;

    free(condWeights);
  }

  PutRNGstate();

  free(newPart);

  end = clock();
  *timing = ((double) (end - start)) / CLOCKS_PER_SEC;

  return;
}

double computeWeightForOneSetBR(int *idxSet, int *nCond, int *partition, double *cov,
				double *sigma2, double *covChol, double *ham, double *mean1,
				double *y){

  /* This function computes the weights of a single set of a given
     partition. Hence the weight of a partition is the product of
     these single weights. (Brown-Resnick case) */
  double weight = 1;
  int r = 0;

  for (int i=0;i<*nCond;i++)
    r += (partition[i] == *idxSet);

  int nr = *nCond - r,
    *tau = malloc(r * sizeof(int)),
    *taubar = malloc(nr * sizeof(int));

  gettau(nCond, partition, idxSet, tau);

  if (r < *nCond) {
    gettaubar(nCond, partition, idxSet, taubar);

    // Build the matrices J and Jtilde
    double *J = malloc(*nCond * nr * sizeof(double)),
      *Jtilde = malloc(*nCond * *nCond * sizeof(double));
    buildJ(tau, nCond, &r, J);
    buildJtilde(tau, nCond, &r, Jtilde);

    /* Define the parameters of the log-normal distribution, i.e.,
       the Cholesky decomposition of the inverse of the covariance
       matrix and the mean mu */
    double *iBchol = malloc(nr * nr * sizeof(double)),
      *mu = malloc(nr * sizeof(double));
    getParametersBR(J, Jtilde, nCond, &nr, covChol, ham, mean1, y, iBchol,
		    mu);

    /* Compute the proba values --- see the note */
    double prob = 0;
    computeprobaBR(iBchol, mu, y, nCond, &r, nCond, taubar, &prob);
    weight = prob;

    free(J); free(Jtilde); free(iBchol); free(mu);
  }

  /* Get the sub-matrix of cov */
  double *covjchol = malloc(r * r * sizeof(double));
  getSubMatrix(cov, nCond, &r, tau, &r, tau, covjchol);

  /* Get the Cholesky decomposition of this submatrix */
  int info = 0;
  F77_CALL(dpotrf)("U", &r, covjchol, &r, &info FCONE);

  if (info != 0)
    error("4. error code %d from Lapack routine '%s'", info, "dpotrf");

  //Get the sub-vector of sigma2 and of y
  double *sigma2j = malloc(r * sizeof(double)),
    *yj = malloc(r * sizeof(double));
  for (int k=0;k<r;k++){
    sigma2j[k] = sigma2[tau[k]];
    yj[k] = y[tau[k]];
  }

  /* Compute the f values --- see the note */
  double f = 0;
  getfvaluesBR(yj, sigma2j, covjchol, &r, &f);
  // f values are OK
  weight *= exp(f);

  free(tau); free(taubar); free(sigma2j); free(yj); free(covjchol);

  return weight;
}

void totoBR(int *nSim, int *n, double *coord, double *ans, double *range,
	  double *smooth){

  GetRNGstate();
  int oneInt = 1, blockSize = 500;
  double *vario = malloc(*n * sizeof(double)),
    *prop = malloc(*n * sizeof(double)),
    *covChol = malloc(*n * *n * sizeof(double)),
    irange = 1 / *range;

  // Build the covariance matrix
  for (int i=0;i<*n;i++){
    for (int j=0;j<=i;j++){
      covChol[i + *n * j] = covChol[j + *n * i] = R_pow(fabs(coord[i]) * irange, *smooth) +
	R_pow(fabs(coord[j]) * irange, *smooth) -
	R_pow(fabs(coord[i] - coord[j]) * irange, *smooth);
    }
  }

  // Compute the Cholesky decomposition
  int info = 0;
  F77_CALL(dpotrf)("U", n, covChol, n, &info FCONE);

  if (info != 0)
    error("Error code %d in Lapack routine '%s'", info, "dpotrf");

  // Compute the variogram
  for (int i=0;i<*n;i++)
    vario[i] = R_pow(fabs(coord[i]) * irange, *smooth);

  for (int i=0;i<*nSim * *n;i++)
    ans[i] = -1e6;

  for (int i=0;i<*nSim;i++){
    double poisson = 0;
    int nBlock = blockSize;

    while (nBlock){
      for (int j=0;j<*n;j++)
	prop[j] = norm_rand();

      F77_CALL(dtrmv)("U", "T", "N", n, covChol, n, prop, &oneInt
	   FCONE FCONE FCONE);

      poisson += exp_rand();

      for (int j=0;j<*n;j++)
	prop[j] = prop[j] - vario[j] - log(poisson);

      for (int j=0; j<*n;j++)
	ans[i + j * *nSim] = fmax2(ans[i + j * *nSim], prop[j]);

      nBlock--;
    }
  }

  PutRNGstate();

  free(vario);free(prop);free(covChol);

  for (int i=0;i<(*n * *nSim);i++)
    ans[i] = exp(ans[i]);

  return;
}


void totoBR2(int *nSim, int *n, double *coord, double *ans, double *range,
	   double *smooth, double *xlim){

  GetRNGstate();
  int oneInt = 1, blockSize = 100, nPoissSup = 200;
  double varioShift[*n], prop[*n], inflate = *range * R_pow(2, 1 / *smooth),
    bounds[2], covShift[*n * *n], irange = 1 / *range;

  bounds[0] = xlim[0] - inflate;
  bounds[1] = xlim[1] + inflate;


  for (int i=0;i<*nSim * *n;i++)
    ans[i] = -1e6;

  for (int i=0;i<*nSim;i++){

    for (int nSup=0;nSup<nPoissSup;nSup++){
      double poisson = 0;
      int nBlock = blockSize;

      while (nBlock){
	double xShift = bounds[0] + (bounds[1] - bounds[0]) * unif_rand();

	for (int j=0;j<*n;j++)
	  varioShift[j] = R_pow(fabs(coord[j] - xShift) * irange, *smooth);

	for (int j=0;j<*n;j++){
	  for (int k=0;k<=j;k++){
	    covShift[j + *n * k] = 0;
	    covShift[k + *n * j] = varioShift[j] + varioShift[k] -
	      R_pow(fabs(coord[j] - coord[k]) * irange, *smooth);
	  }
	}

	// Cholesky decomposition
	int info = 0;
	F77_CALL(dpotrf)("U", n, covShift, n, &info FCONE);

	if (info != 0)
	  error("Error code %d from Lapack routine '%s'", info, "dpotrf");

	for (int j=0;j<*n;j++)
	  prop[j] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", n, covShift, n, prop, &oneInt
	     FCONE FCONE FCONE);

	poisson += nPoissSup * exp_rand();

	for (int j=0;j<*n;j++)
	  prop[j] = prop[j] - varioShift[j] - log(poisson);

	for (int j=0; j<*n;j++)
	  ans[i + j * *nSim] = fmax2(ans[i + j * *nSim], prop[j]);

	nBlock--;
      }
    }

    for (int j=0;j<*n;j++)
      ans[i + j * *nSim] = exp(ans[i + j * *nSim]);

  }

  PutRNGstate();

  return;
}


void getStartingPartitionBR(int *nSim, int *n, double *coord, double *range, double *smooth,
			    int *startPart){

  GetRNGstate();
  int oneInt = 1, blockSize = 500, dummy[*n];
  double vario[*n], prop[*n], covChol[*n * *n], x[*n], irange = 1 / *range;

  // Build the covariance matrix
  for (int i=0;i<*n;i++){
    for (int j=0;j<=i;j++){
      covChol[i + *n * j] = 0;
      covChol[j + *n * i] = R_pow(fabs(coord[i]) * irange, *smooth) +
	R_pow(fabs(coord[j]) * irange, *smooth) -
	R_pow(fabs(coord[i] - coord[j]) * irange, *smooth);
    }
  }

  // Compute the Cholesky decomposition
  int info = 0;
  F77_CALL(dpotrf)("U", n, covChol, n, &info FCONE);

  if (info != 0)
    error("Error code %d in Lapack routine '%s'", info, "dpotrf");

  // Compute the variogram
  for (int i=0;i<*n;i++)
    vario[i] = R_pow(fabs(coord[i]) * irange, *smooth);

  for (int i=0;i<*nSim;i++){
    double poisson = 0;
    int nBlock = blockSize;

    for (int j=0;j<*n;j++){
      x[j] = -1e6;
      dummy[j] = -1;
    }

    int partSet = 0;
    while (nBlock){
      for (int j=0;j<*n;j++)
	prop[j] = norm_rand();

      F77_CALL(dtrmv)("U", "T", "N", n, covChol, n, prop, &oneInt
	   FCONE FCONE FCONE);

      poisson += exp_rand();

      for (int j=0;j<*n;j++)
	prop[j] = prop[j] - vario[j] - log(poisson);

      int hasChanged = 0;
      for (int j=0; j<*n;j++){
	if (prop[j] > x[j]){
	  dummy[j] = partSet;
	  hasChanged = 1;
	}

	x[j] = fmax2(x[j], prop[j]);
      }

      if (hasChanged){
	partSet++;
	convert2rightformat(dummy, n, &partSet);
	partSet = getPartSize(dummy, n);
      }

      nBlock--;
    }

    for (int j=0;j<*n;j++)
      startPart[i * *n + j] = dummy[j];
  }

  return;
}



/* The code below has to be updated for the extremal-t model */

void computeWeightsExtt(int *nCond, double *y, int *nPart, int *allPart,
			int *allSize, double *cov, double *nu, double *weights){
  /* This function computes the distribution function of the random
     partition given the conditioning values, i.e., Pr[\theta = \tau |
     Z(x) = z]. (Extremal-t case) */

  int *currentPart = malloc(*nCond * sizeof(int));

  // Set the weights to one
  for (int i=0;i<*nPart;i++)
    weights[i] = 1;

  for (int i=0;i<*nPart;i++){
    // Loop over all the partitions
    memcpy(currentPart, allPart + i * *nCond, *nCond * sizeof(int));

    for (int j=0;j<allSize[i];j++){
      //Loop over all sets of the current partition

      int r = 0;
      for (int k=0;k<*nCond;k++)
	r += (currentPart[k] == j);

      int nr = *nCond - r,
	*tau = malloc(r * sizeof(int)),
	*taubar = malloc(nr * sizeof(int));

      gettau(nCond, currentPart, &j, tau);

      if (r < *nCond) {
	gettaubar(nCond, currentPart, &j, taubar);

	/* Define the parameters of the distribution of the extremal
	   function. For the Schlather case, it is a multivariate
	   Student distribution up to a deterministic
	   transformation. Hence below we compute the degree of
	   freedom (DoF), the Cholesky decomposition of the scale
	   matrix (scaleMat) and the location parameter (mu) */
	double DoF = 0,
	  *scaleMat = malloc(nr * nr * sizeof(double)),
	  *mu = malloc(nr * sizeof(double));
	getParametersExtt(tau, taubar, &r, &nr, cov, y, nu, &DoF, mu, scaleMat);

	/* Compute the proba values, i.e., the probabilities that one
	   extremal function is smaller than the hitting bounds at the
	   remaining conditioning locations. */
	double prob = 0;
	computeprobaExtt(nu, &DoF, mu, scaleMat, y, &nr, taubar, &prob);
	weights[i] *= prob;

	free(scaleMat); free(mu);
      }

      /* Compute the intensity function values, i.e.,
	 \lambda_{x_{\tau_j}}(z_{\tau_j}). */
      double f = 0;
      getfvaluesExtt(y, nCond, &r, tau, cov, nu, &f);

      weights[i] *= exp(f);
      free(tau); free(taubar);
    }
  }

  /* Compute the normalizing constant and normalize the
     distribution */
  double isumWeights = 0;
  for (int i=0;i<*nPart;i++)
    isumWeights += weights[i];

  isumWeights = 1 / isumWeights;

  for (int i=0;i<*nPart;i++)
    weights[i] *= isumWeights;

  free(currentPart);
  return;
}

void computeprobaExtt(double *nu, double *DoF, double *mu, double *scaleMat, double *y,
		      int *ntaubar, int *taubar, double *prob){

  *prob = 0;

  double *ytaubar = malloc(*ntaubar * sizeof(double));
  for (int i=0;i<*ntaubar;i++)
    ytaubar[i] = R_pow(y[taubar[i]], 1 / *nu);

  if (*ntaubar == 1)
    *prob = pt(*ytaubar - *mu, *DoF, 1, 0);

  else {
    /* Computing the proba */
    double err = 0;
    int nMc = 0;

    pmvt(ytaubar, ntaubar, DoF, mu, scaleMat, prob, &err, &nMc);
  }

  free(ytaubar);

  return;
}

void getfvaluesExtt(double *y, int *n, int *ntau, int *tau, double *cov,
		    double *nu, double *f){

  /* This function computes the logarithm of the intensity
     function, i.e., \lambda_{x_\tau_j}(z_{\tau_j}). */

  /* Get the sub-matrix of cov */
  const double logCnu = (1 - 0.5 * *nu) * M_LN2 + M_LN_SQRT_PI -
    lgammafn(0.5 * (*nu + 1));
  double *covjchol = malloc(*ntau * *ntau * sizeof(double));
  getSubMatrix(cov, n, ntau, tau, ntau, tau, covjchol);

  /* Get the Cholesky decomposition of this submatrix */
  int info = 0;
  F77_CALL(dpotrf)("U", ntau, covjchol, ntau, &info FCONE);

  if (info != 0)
    error("4. error code %d from Lapack routine '%s'", info, "dpotrf");

  //Get the sub-vector of y
  double *yj = malloc(*ntau * sizeof(double));
  for (int k=0;k<*ntau;k++)
    yj[k] = y[tau[k]];

  /* Compute the log of the determinant of cov */
  int oneInt = 1;
  double det = 0;
  for (int i=0;i<*ntau;i++)
    det += log(covjchol[i * (1 + *ntau)]);

  det *= 2;

  /* Compute covjchol^{-T} %*% y */
  double *covjCholMean = malloc(*ntau * sizeof(double));
  for (int i=0;i<*ntau;i++)
    covjCholMean[i] = R_pow(yj[i], 1 / *nu);

  F77_CALL(dtrsv)("U", "T", "N", ntau, covjchol, ntau, covjCholMean, &oneInt
       FCONE FCONE FCONE);

  double mahal = 0;
for (int i=0;i<*ntau;i++){
    mahal += covjCholMean[i] * covjCholMean[i];
    *f += (1 / *nu -  1) * log(yj[i]);
 }

*f += logCnu + (1 - *ntau) * log(*nu) + (2 - *nu) * M_LN2 - *ntau * M_LN_SQRT_PI -
  0.5 * det - 0.5 * (*ntau + *nu) * log(mahal) + lgammafn(0.5 * (*ntau + *nu));

  free(covjCholMean);
  return;
}

void getParametersExtt(int *tau, int *taubar, int *ntau, int *ntaubar, double *cov,
		       double *y, double *nu, double *DoF, double *mu, double *scaleMat){

  /* This function computes the parameters of distribution of the
     extremal functions. For the Extremal-t case, it is a multivariate
     Student distribution with degree of freedom (DoF), location
     parameter (mu) and scale matrix (scaleMat) (but don't forget the
     mapping see the paper). */

  int oneInt = 1;

  // First get the submatrices of Sigma and the subvector of y
  int dim = *ntau + *ntaubar;

  double *covCholx = malloc(*ntau * *ntau * sizeof(double));
  getSubMatrix(cov, &dim, ntau, tau, ntau, tau, covCholx);

  int info = 0;
  F77_CALL(dpotrf)("U", ntau, covCholx, ntau, &info FCONE);
  if (info != 0)
    error("0. error code %d from Lapack routine '%s'", info, "dpotrf");

  double *covs = malloc(*ntaubar * *ntaubar * sizeof(double));
  getSubMatrix(cov, &dim, ntaubar, taubar, ntaubar, taubar, covs);

  double *covsx = malloc(*ntaubar * *ntau * sizeof(double));
  getSubMatrix(cov, &dim, ntaubar, taubar, ntau, tau, covsx);

  double *yx = malloc(*ntau * sizeof(double));
  for (int i=0;i<*ntau;i++)
    yx[i] = y[tau[i]];

  /* Get the degree of freedom parameter */
  *DoF = *ntau + *nu;

  /* Compute the location parameter:

     It is given by mu = covsx %*% covx^(-1) %*% yx, so in terms of
     the Cholesky decomposition we will compute it as follows

     mu = covsx %*% (covCholx^T %*% covCholx)^(-1) %*% (yx)^(1/nu)
        = covsx %*% covCholx^(-1) %*% covCholx^(-T) %*% (yx)^(1/nu),

     i.e.,
       a) solve the (triangular) system X %*% covCholx = covsx -->> dummy1
       b) solve the (triangular) system covCholx^T %*% Y = (yx)^(1/nu)  -->> dummy2
       c) Then mu = X %*% Y.

     Note that Y is a vector and X is a matrix!

  */

  double *dummy1 = malloc(*ntaubar * *ntau * sizeof(double));
  memcpy(dummy1, covsx, *ntaubar * *ntau * sizeof(double));

  double alpha = 1;
  F77_CALL(dtrsm)("R", "U", "N", "N", ntaubar, ntau, &alpha, covCholx, ntau,
		  dummy1, ntaubar
		  FCONE FCONE FCONE FCONE);

  double *dummy2 = malloc(*ntau * sizeof(double));
  for (int i=0;i<*ntau;i++)
    dummy2[i] = R_pow(y[i], 1 / *nu);
  F77_CALL(dtrsv)("U", "T", "N", ntau, covCholx, ntau, dummy2, &oneInt
		  FCONE FCONE FCONE);

  double beta = 0;
  F77_CALL(dgemv)("N", ntaubar, ntau, &alpha, dummy1, ntaubar, dummy2, &oneInt,
		  &beta, mu, &oneInt FCONE);

  /* Compute the scale matrix:

     It is given by

     scaleMat = mahal / (k+*nu) * (covs - covsx %*% covx^(-1) %*% covsx^T),

     so in terms of the Cholesky decomposition we will compute it as
     follows

     dummy3 = covsx %*% covx^(-1) %*% covsx^T
            = covsx %*% (covCholx^T %*% covCholx)^(-1) %*% covsx^T
	    = covsx %*% covCholx^(-1) %*% covCholx^(-T) %*% covsx^T
	    = covsx %*% covCholx^(-1) %*% (covsx %*% covCholx^(-1))^T
	    = dummy1 %*% dummy1^T.

     dummy3 = covs - dummy3;

     mahal = {(yx)^(1/nu)}^T %*% covx^(-1) %*% (yx)^(1/nu)
           = {(yx)^(1/nu)}^T %*% (covCholx^T %*% covCholx)^(-1) %*% (yx)^(1/nu)
	   = (covCholx^(-T) %*% (yx)^(1/nu))^T %*% covCholx^(-T) %*% (yx)^(1/nu)
	   = dummy2^T %*% dummy2

     scaleMat = mahal / (k + nu) * dummy3 = mahal / DoF * dummy3.

  */

  double mahal = 0;
  for (int i=0;i<*ntau;i++)
    mahal += dummy2[i] * dummy2[i];

  mahal /= *DoF;

  double minusMahal = - mahal;

  memcpy(scaleMat, covs, *ntaubar * *ntaubar * sizeof(double));
  F77_CALL(dsyrk)("U", "N", ntaubar, ntau, &minusMahal, dummy1, ntaubar,
		  &mahal, scaleMat, ntaubar
		  FCONE FCONE);// Watch out only the upper part is stored.

  // Fill the lower triangular part
  for (int i=0;i<*ntaubar;i++)
    for (int j=i;j<*ntaubar;j++)
      scaleMat[j + i * *ntaubar] = scaleMat[i + *ntaubar * j];

  free(covCholx); free(covs); free(covsx); free(yx); free(dummy1); free(dummy2);
  return;
}

void condsimextt(int *nsim, int *n, int *nCond, int *allPart, double *nu,
		 double *cov, double *y, double *sim, double *subextfct,
		 double *extfct, double *timings){

  clock_t start, end;
  int oneInt = 1,
    *currentPart = malloc(*nCond * sizeof(int));
  double *x = malloc(*n * sizeof(double)),
    *z = malloc(*n * sizeof(double)),
    *prop = malloc(*n * sizeof(double));

  /* Get the Cholesky decomposition of the covariance matrix */
  int info = 0;
  double *covChol = malloc(*n * *n * sizeof(double));
  memcpy(covChol, cov, *n * *n * sizeof(double));
  F77_CALL(dpotrf)("U", n, covChol, n, &info FCONE);

  if (info != 0)
    error("Error code %d from Lapack routine '%s'", info, "dpotrf");

  GetRNGstate();

  for (int i=0;i<*nsim;i++){

    /* Select the partition */
    memcpy(currentPart, allPart + i * *nCond, *nCond * sizeof(int));
    int size = getPartSize(currentPart, nCond);

    start = clock();
    for (int j=0;j<*n;j++)
      z[j] = -1e6;

    for (int j=0;j<size;j++){

      int r = 0;
      for (int k=0;k<*nCond;k++)
	r += (currentPart[k] == j);

      int nr = *n - r,
	*tau = malloc(r * sizeof(int));

      gettau(nCond, currentPart, &j, tau);

      /* Watch out taubar typically contains some conditioning
	 locations and the locations where we want to get conditional
	 realizations!!! */
      int *taubar = malloc(nr * sizeof(int));
      for (int k=0;k<(*n - *nCond);k++)
	taubar[*nCond - r + k] = *nCond + k;

      if (r < *nCond)
	gettaubar(nCond, currentPart, &j, taubar);

      double DoF = 0,
	*scaleMat = malloc(nr * nr * sizeof(double)),
	*mu = malloc(nr * sizeof(double));
      getParametersExtt(tau, taubar, &r, &nr, cov, y, nu, &DoF, mu, scaleMat);

      // Compute the Cholesky decomposition of scaleMat
      double *scaleMatChol = malloc(nr * nr * sizeof(double));
      memcpy(scaleMatChol, scaleMat, nr * nr * sizeof(double));

      int info = 0;
      F77_CALL(dpotrf)("U", &nr, scaleMatChol, &nr, &info FCONE);

      if (info != 0)
	error("2. error code %d from Lapack routine '%s'", info, "dpotrf");

      /* Simulate student processes that satisfy the upper bound
	 constraints from iBchol */
      for (int k=0;k<r;k++)
	z[tau[k]] = y[tau[k]];

      double *eps = malloc(nr * sizeof(double));
      int flag = 1;

      while (flag){
	flag = *nCond - r;

	double scaleFactor = sqrt(DoF / rchisq(DoF));

	for (int k=0;k<nr;k++)
	  eps[k] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", &nr, scaleMatChol, &nr, eps, &oneInt
			FCONE FCONE FCONE);

	for (int k=0;k<nr;k++)
	  eps[k] = R_pow(fmax2(0, mu[k] + scaleFactor * eps[k]), *nu);

	for (int k=0;k<(*nCond - r);k++)
	  flag -= (eps[k] <= y[taubar[k]]);
      }

      for (int k=0;k<nr;k++)
	z[taubar[k]] = fmax2(z[taubar[k]], eps[k]);

      free(scaleMatChol); free(tau); free(taubar); free(mu); free(scaleMat);
      free(eps);
    }

    end = clock();
    timings[0] = ((double) (end - start)) / CLOCKS_PER_SEC;

    /* Simulate extremal-t processes that satisfy the upper bound
       constraints */
    start = clock();
    const double normCst = M_SQRT_PI * R_pow(2, 1 - 0.5 * *nu) /
      gammafn(0.5 * (*nu + 1));
    double poisson = 0;

    for (int j=0;j<*n;j++)
      x[j] = -1e6;

    int nKO = *n;
    while (nKO){

      poisson += exp_rand();
      double ipoisson = normCst / poisson,
	thresh = 3.5 * normCst * ipoisson;

      for (int j=0;j<*n;j++)
	prop[j] = norm_rand();

      F77_CALL(dtrmv)("U", "T", "N", n, covChol, n, prop, &oneInt
		      FCONE FCONE FCONE);

      for (int j=0;j<*n;j++)
	prop[j] = ipoisson * R_pow(fmax(0, prop[j]), *nu);

      int flag = *nCond;
      for (int j=0;j<*nCond;j++)
	flag -= (prop[j] <= y[j]);

      nKO = *n;
      if (flag == 0){
	for (int j=0;j<*n;j++)
	  x[j] = fmax2(x[j], prop[j]);

	//nKO--;
	for (int j=0;j<*n;j++)
	  nKO -= (thresh <= x[j]);
      }
    }

    end = clock();
    timings[1] = ((double) (end - start)) / CLOCKS_PER_SEC;

    /* Get the max between those two simulations */
    for (int j=0;j<*n;j++){
      sim[i * *n + j] = fmax2(x[j], z[j]);
      subextfct[i * *n + j] = x[j];
      extfct[i * *n + j] = z[j];
    }
  }

  PutRNGstate();

  free(currentPart); free(x); free(z); free(prop);
  return;
}

void gibbsForPartExtt(int *nchain, int *nthin, int *burnin, int *nCond,
		      int *currentPart, double *nu, double *cov, double *y, int *chain,
		      double *timings){

  /* This function generates a Markov chain using a Gibbs sampler
     whose target distribution is \Pr[\theata = \tau | Z(x) = z].  The
     update is performed by picking up randomly one conditional
     location, compute all the weights related to moving this location
     to another set (or even a new one), and sample the next state of
     the chain from the conditional distribution. */

  clock_t start, end;
  start = clock();
  int *newPart = malloc(*nCond * sizeof(int)),
    size = getPartSize(currentPart, nCond),
    iterThin = 0,
    iter = 0;

  GetRNGstate();

  while (iterThin < *nchain){

    void R_CheckUserInterrupt(void);

    /* Pick up one location randomly */
    int idxSite = *nCond * unif_rand(), currentSet = currentPart[idxSite];
    memcpy(newPart, currentPart, *nCond * sizeof(int));

    int r = 0;
    for (int j=0;j<*nCond;j++)
      r += (currentPart[j] == currentSet);

    int ncondWeights = size + (r > 1);
    /* Rmk: If r = 1, then we cannot have newSet = size since this
       partition has already been considered when newSet = oldSet */

    double *condWeights = malloc(ncondWeights * sizeof(double));

    for (int newSet=0;newSet<ncondWeights;newSet++){
      if (newSet == currentSet){
	/* This means that the site is moved to the same partition set
	   --> the partition is left unchanged */
	condWeights[newSet] = 1;
	continue;
      }

      newPart[idxSite] = newSet;

      condWeights[newSet] = computeWeightForOneSetExtt(&newSet, nCond, newPart, nu, cov, y) /
	computeWeightForOneSetExtt(&currentSet, nCond, currentPart, nu, cov, y);

      if (r > 1)
	condWeights[newSet] *= computeWeightForOneSetExtt(&currentSet, nCond, newPart, nu, cov, y);

      int r2 = 0;
      for (int k=0;k<*nCond;k++)
	r2 += (currentPart[k] == newSet);

      if (r2 > 0)
	condWeights[newSet] /= computeWeightForOneSetExtt(&newSet, nCond, currentPart, nu, cov, y);
    }

    /* Compute the normalizing constant and normalize the
       distribution. */
    double isumCondWeights = 0;
    for (int j=0;j<ncondWeights;j++)
      isumCondWeights += condWeights[j];

    isumCondWeights = 1 / isumCondWeights;

    for (int j=0;j<ncondWeights;j++)
      condWeights[j] *= isumCondWeights;

    /* Sample from the full conditional distribution. */
    double u = unif_rand();
    int newSet = -1;
    {
      int flag = 1;
      double cumWeight = 0;
      while (flag){
	newSet++;
	cumWeight += condWeights[newSet];
	flag = (u >= cumWeight);
      }
    }

    /* Update the Markov chain. */
    if (newSet != currentSet){
      currentPart[idxSite] = newSet;

      if (r == 1)
	size--;

      else
	size += (newSet == size);

      convert2rightformat(currentPart, nCond, &size);
    }

    if ((iter > *burnin) & ((iter % *nthin) == 0)){
      memcpy(chain + iterThin * *nCond, currentPart, *nCond * sizeof(int));
      iterThin++;
    }

    iter++;
    free(condWeights);
  }

  PutRNGstate();

  free(newPart);

  end = clock();
  *timings = ((double) (end - start)) / CLOCKS_PER_SEC;

  return;
}

double computeWeightForOneSetExtt(int *idxSet, int *nCond, int *partition,
				  double *nu, double *cov, double *y){

  /* This function computes the weights of a single set of a given
     partition. Hence the weight of a partition is the product of
     these single weights. (Schlather case) */
  double weight = 1;
  int r = 0;

  for (int i=0;i<*nCond;i++)
    r += (partition[i] == *idxSet);

  int nr = *nCond - r,
    *tau = malloc(r * sizeof(int)),
    *taubar = malloc(nr * sizeof(int));

  gettau(nCond, partition, idxSet, tau);

  if (r < *nCond) {
    gettaubar(nCond, partition, idxSet, taubar);

    /* Define the parameters of the Student distribution, i.e., the
       degree of freedom (DoF), the Cholesky decomposition of the
       (scaleMat) scale matrix and the location parameter (mu) */
    double DoF = 0,
      *scaleMat = malloc(nr * nr * sizeof(double)),
      *mu = malloc(nr * sizeof(double));
    getParametersExtt(tau, taubar, &r, &nr, cov, y, nu, &DoF, mu, scaleMat);

    /* Compute the proba values --- see the note */
    double prob = 0;
    computeprobaExtt(nu, &DoF, mu, scaleMat, y, &nr, taubar, &prob);
    weight = prob;

    free(scaleMat); free(mu);
  }

  /* Compute the f values --- see the note */
  double f = 0;
  getfvaluesExtt(y, nCond, &r, tau, cov, nu, &f);
  // f values are OK

  weight *= exp(f);

  free(tau); free(taubar);
  return weight;
}

void totoExtt(int *nsim, int *n, double *nu, double *covChol, double *ans){

  int oneInt = 1;
  const double normCst = M_SQRT_PI * R_pow(2, 1 - 0.5 * *nu) /
    gammafn(0.5 * (*nu + 1));

  double *x = malloc(*n * sizeof(double)),
    *prop = malloc(*n * sizeof(double));

  GetRNGstate();
  for (int i=0;i<*nsim;i++){
    for (int j=0;j<*n;j++)
      x[j] = 0;

    double poisson = 0;
    int nKO = *n;

    while(nKO){
      poisson += exp_rand();

      double ipoisson = 1 / poisson,
	thresh = 3.5 * normCst * ipoisson;

      for (int j=0;j<*n;j++)
	prop[j] = norm_rand();

      F77_CALL(dtrmv)("U", "T", "N", n, covChol, n, prop, &oneInt
		      FCONE FCONE FCONE);

      for (int j=0;j<*n;j++)
	prop[j] = normCst * ipoisson * R_pow(fmax2(0, prop[j]), *nu);

      for (int j=0;j<*n;j++)
	x[j] = fmax2(x[j], prop[j]);

      nKO = *n;
      for (int j=0;j<*n;j++)
	nKO -= (thresh <= x[j]);
    }

    for (int j=0;j<*n;j++)
      ans[i + j * *nsim] = x[j];
  }
  PutRNGstate();
  free(x); free(prop);

  return;
}

void getStartingPartitionExtt(int *nsim, int *n, double *nu, double *covChol,
			      int *startPart){

  /* This function simulates a bunch of max-stable realizations and
     define a suitable starting values for the Gibbs sampler. */
  int oneInt = 1,
    *dummy = malloc(*n * sizeof(int));
  double *x = malloc(*n * sizeof(double)),
    *prop = malloc(*n * sizeof(double));

  for (int i=0;i<*nsim;i++){
    for (int j=0;j<*n;j++){
      x[j] = 0;
      dummy[j] = -1;
    }

    double poisson = 0;
    int nKO = *n;

    GetRNGstate();
    int partSet = 0;
    while(nKO){
      poisson += exp_rand();

      double ipoisson = 1 / poisson,
	thresh = 3.5 * ipoisson;

      for (int j=0;j<*n;j++)
	prop[j] = norm_rand();

      F77_CALL(dtrmv)("U", "T", "N", n, covChol, n, prop, &oneInt
		      FCONE FCONE FCONE);

      for (int j=0;j<*n;j++)
	prop[j] = ipoisson * R_pow(fmax2(0, prop[j]), *nu);

      int hasChanged = 0;
      for (int j=0;j<*n;j++){
	if (prop[j] > x[j]){
	  dummy[j] = partSet;
	  hasChanged = 1;
	}

	x[j] = fmax2(x[j], prop[j]);
      }

      nKO = *n;
      for (int j=0;j<*n;j++)
	nKO -= (thresh <= x[j]);

      if (hasChanged){
	partSet++;
	convert2rightformat(dummy, n, &partSet);
	partSet = getPartSize(dummy, n);
      }
    }

    for (int j=0;j<*n;j++)
      startPart[i * *n + j] = dummy[j];
  }

  PutRNGstate();

  free(x); free(prop);

  return;
}
