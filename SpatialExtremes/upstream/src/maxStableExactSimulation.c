#include "header.h"

void rbrownexact(double *coord, int *nObs, int *nSite, int *dim,
		 int *grid, double *range, double *smooth,
		 double *ans){

  /*
    This function generates random fields for the Brown-Resnick model
    using the exact procedure of Dombry et al. (2106) "Exact
    simulation of max-stable processes" Biometrika

    coord: the coordinates of the locations
    nObs: the number of observations to be generated
    nSite: the number of locations
    dim: the random field is generated in R^dim
    grid: Does coord specifies a grid?
    range: the range parameter
    smooth: the smooth parameter
    ans: the generated random field
  */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1;
  double zero = 0, one = 1, irange = 1 / *range;
  int covmod = 6;//i.e, fractional Brownian motion

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double)),
    *vario = malloc(neffSite * sizeof(double)),
    *shiftedCoord = malloc(*nSite * *dim * sizeof(double)),
    *orig = malloc(*dim * sizeof(double)),
    *poisson = malloc(*nObs * sizeof(double));

  buildcovmat(nSite, grid, &covmod, coord, dim, &zero, &one, range,
	      smooth, covmat);

  /* Compute the Cholesky decomposition of the covariance matrix once for all */
  int info = 0;
  F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info FCONE);

  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dpotrf");

  GetRNGstate();
  for (int j=0;j<neffSite;j++){
    // Set the origin
    if (*grid){
      int idx1 = j / *nSite, idx2 = j % *nSite;//works only for 2d grid
      orig[0] = coord[idx1];
      orig[1] = coord[*nSite + idx2];
    }
    else {
      for (int d=0;d<*dim;d++)
	orig[d] = coord[j + d * *nSite];
    }

    // Compute the variogram gamma(s - origin)
    for (int l=0; l<*nSite;l++)
      for (int d=0; d<*dim; d++)
	shiftedCoord[d * *nSite + l] = coord[d * *nSite + l] - orig[d];

    distance2orig(shiftedCoord, *nSite, *dim, vario, *grid);

    for (int l=0; l<neffSite; l++)
      vario[l] = R_pow(vario[l] * irange, *smooth);

    for (int i=0; i<*nObs; i++){

      poisson[i] = exp_rand();
      double ipoisson = -log(poisson[i]);

      while (ans[j * lagj + i * lagi] < ipoisson){
	R_CheckUserInterrupt();

	// Generate a proposal extremal function
	for (int l=0;l<neffSite;l++)
	  gp[l] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt

			FCONE FCONE FCONE);

	double dummy = gp[j];
	for (int l=0;l<neffSite;l++)
	  gp[l] -= dummy + vario[l];

	// Update the max-stable realization (if any)
	int valid = 1;
	for (int l=0; l<j; l++){
	  if ((ipoisson + gp[l]) > ans[l * lagj + i * lagi]){
	    valid = 0;
	    break;
	  }
	}

	if (valid == 1)//Valid extremal function --> update \eta(s)
	  for (int l=j;l<neffSite;l++)
	    ans[l * lagj + i * lagi] = fmax2(ans[l * lagj + i * lagi], ipoisson + gp[l]);

	// Update the "norm" of the spectral function (log-scale)
	poisson[i] += exp_rand();
	ipoisson = -log(poisson[i]);
      }
    }
  }

  // Switch to unit Frechet margins
  for (int i=0;i<(*nObs * neffSite);i++)
    ans[i] = exp(ans[i]);

  PutRNGstate();
  free(covmat); free(gp); free(vario); free(shiftedCoord); free(orig); free(poisson);

  return;
}


void rextremaltexact(double *coord, int *nObs, int *nSite, int *dim,
		     int *covmod, int *grid, double *nugget, double *range,
		     double *smooth, double *DoF, double *ans){

  /*
    This function generates random fields for the extremal-t model
    using the exact procedure of Dombry et al. (2106) "Exact
    simulation of max-stable processes" Biometrika

    coord: the coordinates of the locations
    nObs: the number of observations to be generated
    nSite: the number of locations
    dim: the random field is generated in R^dim
    covmod: the covariance model
    grid: Does coord specifies a grid?
    range: the range parameter
    smooth: the smooth parameter
    DoF: the degree of freedom
    ans: the generated random field
  */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1;
  double sill = 1 - *nugget,
    miDoF = - 1.0 / *DoF;

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *scalemat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double)),
    *mu = malloc(neffSite * sizeof(double)),
    *poisson = malloc(*nObs * sizeof(double));

  buildcovmat(nSite, grid, covmod, coord, dim, nugget, &sill, range,
	      smooth, covmat);

  GetRNGstate();
  for (int j=0;j<neffSite;j++){
    // Get the mean vector of the extremal function
    for (int l=0;l<neffSite;l++)
      mu[l] = covmat[l + j * neffSite];

    // Compute the scale matrix of the extremal function
    for (int l=0;l<neffSite;l++)
      for (int m=l;m<neffSite;m++)
	scalemat[l + m * neffSite] = scalemat[m + l * neffSite] =
	  (covmat[l + m * neffSite] - covmat[j + l * neffSite] * covmat[j + m * neffSite]) /
	  (1.0 + *DoF);

    // By construction this matrix is singular since the j-th row is 0
    // so we regularize it
    scalemat[j + j * neffSite] = 1e-12;

    // Compute the Cholesky decomposition of this matrix
    int info = 0;
    F77_CALL(dpotrf)("U", &neffSite, scalemat, &neffSite, &info FCONE);

    if (info != 0)
      error("error code %d from Lapack routine '%s'", info, "dpotrf");

    // Reset it to the degenerate case
    scalemat[j + j * neffSite] = 0;

    for (int i=0; i<*nObs; i++){
      poisson[i] = exp_rand();
      double ipoissonDoF = R_pow(poisson[i], miDoF);

      while (ans[j * lagj + i * lagi] < ipoissonDoF){
	R_CheckUserInterrupt();

	// Generate a proposal extremal function (log-scale)
	for (int l=0;l<neffSite;l++)
	  gp[l] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", &neffSite, scalemat, &neffSite, gp, &oneInt
			FCONE FCONE FCONE);

	double scale = sqrt((1 + *DoF) / rchisq(1 + *DoF));
	for (int l=0;l<neffSite;l++)
	  gp[l] = mu[l] + gp[l] * scale;

	// Update the max-stable realization (if any)
	int valid = 1;
	for (int l=0; l<j; l++){
	  if ((ipoissonDoF * gp[l]) > ans[l * lagj + i * lagi]){
	    valid = 0;
	    break;
	  }
	}

	if (valid == 1)//Valid extremal function --> update \eta(s)
	  for (int l=j;l<neffSite;l++)
	    ans[l * lagj + i * lagi] = fmax2(ans[l * lagj + i * lagi], gp[l] * ipoissonDoF);

	// Update the "norm" of the spectral function (1/DoF scale)
	poisson[i] += exp_rand();
	ipoissonDoF = R_pow(poisson[i], miDoF);
      }
    }
  }

  // Switch to unit Frechet margins
  for (int i=0;i<(*nObs * neffSite);i++)
    ans[i] = R_pow(ans[i], *DoF);

  PutRNGstate();
  free(covmat); free(scalemat); free(gp); free(mu); free(poisson);

  return;
}

void rschlatherexact(double *coord, int *nObs, int *nSite, int *dim,
		     int *covmod, int *grid, double *nugget, double *range,
		     double *smooth, double *ans){

  /*
    This function generates random fields for the Schlather model
    using the exact procedure of Dombry et al. (2106) "Exact
    simulation of max-stable processes" Biometrika

    coord: the coordinates of the locations
    nObs: the number of observations to be generated
    nSite: the number of locations
    dim: the random field is generated in R^dim
    covmod: the covariance model
    grid: Does coord specifies a grid?
    range: the range parameter
    smooth: the smooth parameter
    ans: the generated random field
  */

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
    *scalemat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double)),
    *mu = malloc(neffSite * sizeof(double)),
    *poisson = malloc(*nObs * sizeof(double));

  buildcovmat(nSite, grid, covmod, coord, dim, nugget, &sill, range,
	      smooth, covmat);

  GetRNGstate();
  for (int j=0;j<neffSite;j++){
    // Get the mean vector of the extremal function
    for (int l=0;l<neffSite;l++)
      mu[l] = covmat[l + j * neffSite];

    // Compute the scale matrix of the extremal function
    for (int l=0;l<neffSite;l++)
      for (int m=l;m<neffSite;m++)
	scalemat[l + m * neffSite] = scalemat[m + l * neffSite] =
	  (covmat[l + m * neffSite] - covmat[j + l * neffSite] * covmat[j + m * neffSite]) / 2.0;

    // By construction this matrix is singular since the j-th row is 0
    // so we regularize it
    scalemat[j + j * neffSite] = 1e-12;

    // Compute the Cholesky decomposition of this matrix
    int info = 0;
    F77_CALL(dpotrf)("U", &neffSite, scalemat, &neffSite, &info FCONE);

    if (info != 0)
      error("error code %d from Lapack routine '%s'", info, "dpotrf");

    // Reset it to the degenerate case
    scalemat[j + j * neffSite] = 0;

    for (int i=0; i<*nObs; i++){
      poisson[i] = exp_rand();

      while ((poisson[i] * ans[j * lagj + i * lagi]) < 1){
	R_CheckUserInterrupt();

	// Generate a proposal extremal function (log-scale)
	for (int l=0;l<neffSite;l++)
	  gp[l] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", &neffSite, scalemat, &neffSite, gp, &oneInt
	     FCONE FCONE FCONE);

	double scale = sqrt(2 / rchisq(2));
	for (int l=0;l<neffSite;l++)
	  gp[l] = mu[l] + gp[l] * scale;

	// Update the max-stable realization (if any)
	int valid = 1;
	for (int l=0; l<j; l++){
	  if (gp[l] > (poisson[i] * ans[l * lagj + i * lagi])){
	    valid = 0;
	    break;
	  }
	}

	if (valid == 1)//Valid extremal function --> update \eta(s)
	  for (int l=j;l<neffSite;l++)
	    ans[l * lagj + i * lagi] = fmax2(ans[l * lagj + i * lagi], gp[l] / poisson[i]);

	// Update the "norm" of the spectral function
	poisson[i] += exp_rand();
      }
    }
  }

  PutRNGstate();
  free(covmat); free(scalemat); free(gp); free(mu); free(poisson);

  return;
}

