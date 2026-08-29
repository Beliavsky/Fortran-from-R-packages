#include "header.h"

void buildcovmat(int *nSite, int *grid, int *covmod, double *coord, int *dim,
		 double *nugget, double *sill, double *range,
		 double *smooth, double *covMat){

  int nPairs, effnSite = *nSite, zero = 0;
  const double one = 1, dzero = 0;
  double flag = 0;

  if (*grid)
    effnSite = R_pow_di(effnSite, *dim);

  nPairs = effnSite * (effnSite - 1) / 2;

  double *dist = malloc(nPairs * sizeof(double)),
    *rho = malloc(nPairs * sizeof(double)),
    *coordGrid = malloc(effnSite * *dim * sizeof(double));

  if (*grid){
    //Coord specify a grid
    for (int i = 0; i < *nSite; i++)
      for (int j = 0; j < *nSite; j++){
	coordGrid[j + i * *nSite] = coord[i];
	coordGrid[*nSite * (*nSite + i) + j] = coord[j];
      }

    distance(coordGrid, dim, &effnSite, &zero, dist);
  }

  else
    //Coord don't specify a grid
    distance(coord, dim, nSite, &zero, dist);

  switch (*covmod){
  case 1:
    flag = whittleMatern(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 2:
    flag = cauchy(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 3:
    flag = powerExp(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 4:
    flag = bessel(dist, nPairs, *dim, dzero, one, *range, *smooth, rho);
    break;
  case 6:
    if (*grid)
      flag = fbm(coordGrid, dist, *dim, effnSite, one, *range, *smooth, rho);

    else
      flag = fbm(coord, dist, *dim, effnSite, one, *range, *smooth, rho);

    break;
  }

  if (flag != 0.0)
    error("The covariance parameters seem to be ill-defined. Please check\n");

  //Fill the non-diagonal elements of the covariance matrix
  //#pragma omp parallel for
  for (int currentPair=0;currentPair<nPairs;currentPair++){
    int i = 0, j = 0;
    getSiteIndex(currentPair, effnSite, &i, &j);
    covMat[effnSite * i + j] = covMat[effnSite * j + i] = *sill * rho[currentPair];
  }

  //Fill the diagonal elements of the covariance matrix
  if (*covmod == 6){
    //Fractional brownian
    double irange2 = 1 / (*range * *range);

    if (*grid){
      for (int i = 0; i < effnSite;i++){
	covMat[i * (effnSite + 1)] = 0;

	for (int j= 0; j < *dim; j++)
	  covMat[i * (effnSite + 1)] += coordGrid[i + j * effnSite] * coordGrid[i + j * effnSite];

	covMat[i * (effnSite + 1)] = 2 * pow(covMat[i * (effnSite + 1)] * irange2, 0.5 * *smooth);
      }
    }

    else {
      for (int i = 0; i < effnSite; i++){
	covMat[i * (effnSite + 1)] = 0;

	for (int j = 0; j < *dim; j++)
	  covMat[i * (effnSite + 1)] += coord[i + j * effnSite] * coord[i + j * effnSite];

	covMat[i * (effnSite + 1)] = 2 * pow(covMat[i * (effnSite + 1)] * irange2, 0.5 * *smooth);
      }
    }
  }

  else
    for (int i = 0; i < effnSite; i++)
      covMat[i * (effnSite + 1)] = *sill + *nugget;


  free(dist); free(rho); free(coordGrid);
  return;
}

void direct(int *n, int *nSite, int *grid, int *covmod, double *coord, int *dim,
	    double *nugget, double *sill, double *range, double *smooth,
	    double *ans){

  int neffSite = *nSite, lagi = 1, lagj = 1;

  if (*grid){
    neffSite = R_pow_di(neffSite, *dim);
    lagi = neffSite;
  }

  else
    lagj = *n;

  double *covmat = malloc(neffSite * neffSite * sizeof(double));

  buildcovmat(nSite, grid, covmod, coord, dim, nugget, sill, range,
	      smooth, covmat);

  /* Compute the Cholesky decomposition of the covariance matrix */
  int info = 0;
  F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info FCONE);

  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dpotrf");

  /* Simulation part */
  GetRNGstate();

  for (int i=0;i<*n;i++){
    for (int j=0;j<neffSite;j++)
      ans[j * lagj + i * lagi] = norm_rand();

    F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite,
		    ans + i * lagi, &lagj FCONE FCONE FCONE);
  }

  PutRNGstate();

  free(covmat);
  return;
}


