#include "header.h"

void rextremalttbm(double *coord, int *nObs, int *nSite, int *dim,
		   int *covmod, int *grid, double *nugget, double *range,
		   double *smooth, double *DoF, double *uBound, int *nlines,
		   double *ans){
  /* This function generates random fields from the Extremal-t model

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
    nSite: the number of locations
       dim: the random field is generated in R^dim
    covmod: the covariance model
      grid: Does coord specifies a grid?
      nugget: the nugget parameter
     range: the range parameter
    smooth: the smooth parameter
       DoF: the degree of freedom
 blockSize: simulated field is the maximum over blockSize ind. replicates
    nlines: the number of lines used for the TBM algo
       ans: the generated random field */

  int i, neffSite, lagi = 1, lagj = 1;
  double sill = 1 - *nugget;
  const double irange = 1 / *range;

  //rescale the coordinates
  for (i=(*nSite * *dim);i--;)
        coord[i] = coord[i] * irange;

  double *lines = malloc(3 * *nlines * sizeof(double));
  
  if ((*covmod == 3) && (*smooth == 2))
    //This is the gaussian case
    *covmod = 5;

  //Generate lines
  vandercorput(nlines, lines);

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *gp = malloc(neffSite * sizeof(double));

  GetRNGstate();
  
  for (i=*nObs;i--;){
    int nKO = neffSite;
    double poisson = 0;

    while (nKO){

      /* ------- Random rotation of the lines ----------*/
      double u = unif_rand() - 0.5,
	v = unif_rand() - 0.5,
	w = unif_rand() - 0.5,
	angle = runif(0, M_2PI),	
	inorm = 1 / sqrt(u * u + v * v + w * w);
      
      u *= inorm;
      v *= inorm;
      w *= inorm;
      
      rotation(lines, nlines, &u, &v, &w, &angle);
      /* -------------- end of rotation ---------------*/
      
      poisson += exp_rand();
      double ipoisson = 1 / poisson,
	thresh = *uBound * ipoisson;
      
      /* We simulate one realisation of a gaussian random field with
	 the required covariance function */
      for (int j=neffSite;j--;)
	gp[j] = 0;

      tbmcore(nSite, &neffSite, dim, covmod, grid, coord, nugget,
	      &sill, range, smooth, nlines, lines, gp);
      
      nKO = neffSite;
      for (int j=neffSite;j--;){
	double dummy = R_pow(fmax2(0, gp[j]), *DoF) * ipoisson;
	ans[j * lagj + i * lagi] = fmax2(dummy, ans[j * lagj + i * lagi]);
	nKO -= (thresh <= ans[j * lagj + i * lagi]);
      }    
    }
  }
  
  PutRNGstate();
  
  //Lastly we multiply by the normalizing constant
  const double imean = M_SQRT_PI * R_pow(2, -0.5 * (*DoF - 2)) /
    gammafn(0.5 * (*DoF + 1));
  
  for (i=(neffSite * *nObs);i--;)
    ans[i] *= imean;
  
  free(lines); free(gp);
  return;
}

void rextremaltdirect(double *coord, int *nObs, int *nSite, int *dim,
		      int *covmod, int *grid, double *nugget, double *range,
		      double *smooth, double *DoF, double *uBound, double *ans){
  /* This function generates random fields for the Extremal-t model

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
    nSite: the number of locations
       dim: the random field is generated in R^dim
    covmod: the covariance model
      grid: Does coord specifies a grid?
      nugget: the nugget parameter
     range: the range parameter
    smooth: the smooth parameter
       DoF: the degree of freedom
 blockSize: see rextremalttbm.
       ans: the generated random field */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1;
  double sill = 1 - *nugget;

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double));

  buildcovmat(nSite, grid, covmod, coord, dim, nugget, &sill, range,
	      smooth, covmat);
  
  /* Compute the Cholesky decomposition of the covariance matrix */
  int info = 0;
  F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info FCONE);

  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dpotrf");
  
  GetRNGstate();
 
  for (int i=*nObs;i--;){
    double poisson = 0;
    int nKO = neffSite;
      
    while (nKO){
      poisson += exp_rand();

      double ipoisson = 1 / poisson,
	thresh = *uBound * ipoisson;

      
      /* We simulate one realisation of a gaussian random field with
	 the required covariance function */
      for (int j=neffSite;j--;)
	gp[j] = norm_rand();
      
      F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt
		      FCONE FCONE FCONE);
      
      nKO = neffSite;
      for (int j=neffSite;j--;){
	double dummy = R_pow(fmax2(0, gp[j]), *DoF) * ipoisson;
	ans[j * lagj + i * lagi] = fmax2(dummy, ans[j * lagj + i * lagi]);
	nKO -= (thresh <= ans[j * lagj + i * lagi]);
      }
    }
  }

  PutRNGstate();
  //Lastly we multiply by the normalizing constant
  const double imean = M_SQRT_PI * R_pow(2, -0.5 * (*DoF - 2)) /
    gammafn(0.5 * (*DoF + 1));
  for (int i=(neffSite * *nObs);i--;)
    ans[i] *= imean;
  
  free(covmat); free(gp);
  return;
}

void rextremaltcirc(int *nObs, int *ngrid, double *steps, int *dim,
		    int *covmod, double *nugget, double *range,
		    double *smooth, double *DoF, double *uBound, double *ans){
  /* This function generates random fields from the Schlather model

     nObs: the number of observations to be generated
    ngrid: the number of locations along one axis
      dim: the random field is generated in R^dim
   covmod: the covariance model
     nugget: the nugget parameter
    range: the range parameter
   smooth: the smooth parameter
      DoF: the degree of freedom
blockSize: see rextremalttbm
      ans: the generated random field */

  int i, j, k = -1, nbar = R_pow_di(*ngrid, *dim), r, m;
  const double zero = 0;
  double *rho, *irho, sill = 1 - *nugget;
    //Below is a table of highly composite numbers
  int HCN[39] = {1, 2, 4, 6, 12, 24, 36, 48, 60, 120, 180, 240,
		 360, 720, 840, 1260, 1680, 2520, 5040, 7560,
		 10080, 15120, 20160, 25200, 27720, 45360, 50400,
		 55440, 83160, 110880, 166320, 221760, 277200,
		 332640, 498960, 554400, 665280, 720720, 1081080};

    
  /* Find the smallest size m for the circulant embedding matrix */
  {
    int dummy = 2 * (*ngrid - 1);
    do {
      k++;
      m = HCN[k];
    } while (m < dummy);
  }
  
  /* ---------- beginning of the embedding stage ---------- */
  int mbar = m * m, halfM = m / 2, notPosDef = 0;
  do {
    double *dist = (double *)R_alloc(mbar, sizeof(double));

    notPosDef = 0;
    //Computation of the distance
    for (r=mbar;r--;){
      i = r % m;
      j = r / m;
      
      if (i > halfM)
	i -= m;
      
      if (j > halfM)
	j -= m;
      
      dist[r] = hypot(steps[0] * i, steps[1] * j);
    }

    //Computations of the covariances
    rho = (double *)R_alloc(mbar, sizeof(double));
    irho = (double *)R_alloc(mbar, sizeof(double));
    for (i=mbar;i--;)
      irho[i] = 0;

    switch (*covmod){
    case 1:
      whittleMatern(dist, mbar, zero, sill, *range, *smooth, rho);
      break;
    case 2:
      cauchy(dist, mbar, zero, sill, *range, *smooth, rho);
      break;
    case 3:
      powerExp(dist, mbar, zero, sill, *range, *smooth, rho);
      break;
    case 4:
      bessel(dist, mbar, *dim, zero, sill, *range, *smooth, rho);
      break;
    }

    /* Compute the eigen values to check if the circulant embbeding
       matrix is positive definite */

    /* Note : The next lines is only valid for 2d random fields. I
       need to change if there are m_1 \neq m_2 as I suppose that m_1
       = m_2 = m */
    int maxf, maxp, *iwork;
    double *work;

    fft_factor(m, &maxf, &maxp);
    work = (double *)R_alloc(4 * maxf, sizeof(double));
    iwork = (int *)R_alloc(maxp, sizeof(int));
    fft_work(rho, irho, m, m, 1, -1, work, iwork);

    fft_factor(m, &maxf, &maxp);
    work = (double *)R_alloc(4 * maxf, sizeof(double));
    iwork = (int *)R_alloc(maxp, sizeof(int));
    fft_work(rho, irho, 1, m, m, -1, work, iwork);

    //Check if the eigenvalues are all positive
    for (i=mbar;i--;){
      notPosDef |= (rho[i] <= 0) || (fabs(irho[i]) > 0.001);
    }

    if (notPosDef){
      k++;
      m = HCN[k];
      halfM = m / 2;
      mbar = m * m;
    }

    if (k > 30)
      error("Impossible to embbed the covariance matrix");
    
  } while (notPosDef);
  /* --------- end of the embedding stage --------- */

  /* Computation of the square root of the eigenvalues */
  for (i=mbar;i--;){
    rho[i] = sqrt(rho[i]);
    irho[i] = 0;//No imaginary part
  }

  int mdag = m / 2 + 1, mdagbar = mdag * mdag;
  double isqrtMbar = 1 / sqrt(mbar);

  double *a = malloc(mbar * sizeof(double)),
    *ia = malloc(mbar * sizeof(double)),
    *gp = malloc(nbar * sizeof(double));

  GetRNGstate();
  for (int i=*nObs;i--;){
    int nKO = nbar;
    double poisson = 0;

    while (nKO){
      poisson += exp_rand();
      double ipoisson = 1 / poisson,
	thresh = *uBound * ipoisson;
      
      /* We simulate one realisation of a gaussian random field with
	 the required covariance function */
      circcore(rho, a, ia, m, halfM, mdag, mdagbar, *ngrid, nbar, isqrtMbar, *nugget, gp);
      
      nKO = nbar;
      for (int j=nbar;j--;){
	double dummy = R_pow(fmax2(gp[j], 0), *DoF) * ipoisson;
	ans[j + i * nbar] = fmax2(dummy, ans[j + i * nbar]);
	nKO -= (thresh <= ans[j + i * nbar]);
      }
    }
  }
  
  PutRNGstate();
  
  //Lastly we multiply by the normalizing constant
  const double imean = M_SQRT_PI * R_pow(2, -0.5 * (*DoF - 2)) /
    gammafn(0.5 * (*DoF + 1));
  for (i=(nbar * *nObs);i--;)
    ans[i] *= imean;
  
  free(a); free(ia); free(gp);
  return;
}
