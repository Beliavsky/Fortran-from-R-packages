#include "header.h"

void rhitscenbrown(double *coord, int *nObs, int *nSite, int *dim,
		   int *grid, double *range, double *smooth,
		   double *ans, int *ans2){

  /*
    This function generates hitting scenario for the Brown-Resnick model
    using the exact procedure of Dombry et al. (2106) "Exact
    simulation of max-stable processes" Biometrika

    coord: the coordinates of the locations
    nObs: the number of observations to be generated
    nSite: the number of locations
    dim: the random field is generated in R^dim
    grid: Does coord specifies a grid?
    range: the range parameter
    smooth: the smooth parameter
    ans: the generated max-stable process
    ans2: the generated hitting scenario
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

  // Initialize the proposal extremal function index
  int *extFctIdx = malloc(*nObs * sizeof(int));
  for (int i=0;i<*nObs;i++)
    extFctIdx[i] = 0;

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

	F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt FCONE FCONE FCONE);

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

	if (valid == 1){//Valid extremal function --> update \eta(s)
	  extFctIdx[i]++;
	  for (int l=j;l<neffSite;l++){
	    if ((ipoisson + gp[l]) > ans[l * lagj + i * lagi]){
	      ans[l * lagj + i * lagi] = ipoisson + gp[l];
	      ans2[l * lagj + i * lagi] = extFctIdx[i];
	    }
	  }
	}

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
  free(extFctIdx);

  return;
}
