#include "header.h"

void smithstderr(double *data, double *distVec, int *nSite, int *nObs, double *locdsgnmat,
		 int *nloccoeff, double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		 int *nshapecoeff, double *tempdsgnmatloc, int *ntemploccoeff,
		 double *tempdsgnmatscale, int *ntempscalecoeff, double *tempdsgnmatshape,
		 int *ntempshapecoeff, double *loccoeff, double *scalecoeff, double *shapecoeff,
		 double *temploccoeff, double *tempscalecoeff, double *tempshapecoeff,
		 double *cov11, double *cov12, double *cov22, int *fitmarge, int *usetempcov,
		 double *weights, double *hess, double *grad){

  //This is the Smith model. It computes the hessian of the pairwise log-likelihood

  const int nPairs = *nSite * (*nSite - 1) / 2,
    flag = usetempcov[0] + usetempcov[1] + usetempcov[2];
  const double det = *cov11 * *cov22 - *cov12 * *cov12;
  int i, currentPair = -1;
  double  *mahalDist, *locs, *scales, *shapes, *jac, *frech, *trendlocs,
    *trendscales, *trendshapes;

  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  trendlocs = (double *)R_alloc(*nObs, sizeof(double));
  trendscales = (double *)R_alloc(*nObs, sizeof(double));
  trendshapes = (double *)R_alloc(*nObs, sizeof(double));

  /* We have to intialized them to 0 as we don't know if
     dsgnmat2temptrend will be called */
  for (i=*nObs;i--;)
    trendlocs[i] = trendscales[i] = trendshapes[i] = 0;

  //Computing the Mahalanobis distance
  mahalDistFct(distVec, nPairs, cov11, cov12, cov22, mahalDist);

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){

    dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff, shapecoeff,
		  *nSite, *nloccoeff, *nscalecoeff, *nshapecoeff, locs, scales, shapes);

    //Stage 1: Transformation to unit Frechet
    if (flag) {
      dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, temploccoeff,
			tempscalecoeff, tempshapecoeff, *nSite, *nObs, usetempcov, *ntemploccoeff,
			*ntempscalecoeff, *ntempshapecoeff, trendlocs, trendscales, trendshapes);

      gev2frechTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs, trendscales,
		     trendshapes, jac, frech);
    }

    else
      gev2frech(data, *nObs, *nSite, locs, scales, shapes, jac, frech);
  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];

  //Stage 2: Hessian computations
  // a- Covariance matrix part
  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){
      currentPair++;

      int k;
      double imahal = 1 / mahalDist[currentPair], imahalSquare = imahal * imahal;

      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    hess[k * nPairs + currentPair] = hess[(*nObs + k) * nPairs + currentPair] =
	      hess[(2 * *nObs + k) * nPairs + currentPair] = NA_REAL;

	    continue;
	  }

	  double ifrech1 = 1 / frech[k + i * *nObs], ifrech2 = 1 / frech[k + j * *nObs],
	    ifrech1Square = ifrech1 * ifrech1, ifrech2Square = ifrech2 * ifrech2,
	    c1 = log(frech[k + j * *nObs] * ifrech1) * imahal + 0.5 * mahalDist[currentPair],
	    c2 = mahalDist[currentPair] - c1,
	    dnormc1 = dnorm(c1, 0., 1., 0), pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0), pnormc2 = pnorm(c2, 0., 1., 1, 0),
	    //A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	    B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	    dnormc2 * imahal * ifrech2Square,
	    C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	    dnormc1 * imahal * ifrech1Square,
	    D = c2 * dnormc1 * ifrech2 * imahalSquare * ifrech1Square +
	    c1 * dnormc2 * ifrech1 * imahalSquare * ifrech2Square,
	    dAa = - c2 * dnormc1 * ifrech1 * imahal - c1 * dnormc2 * imahal * ifrech2,
	    dBa = (c1 * c1 - 1) * dnormc2 * imahalSquare * ifrech2Square +
	    (1 + c1 * c2 ) * dnormc1 * ifrech1 * ifrech2 * imahalSquare,
	    dCa = (c2 * c2 - 1) * dnormc1 * imahalSquare * ifrech1Square +
	    (1 + c1 * c2) * dnormc2 * ifrech1 * ifrech2 * imahalSquare,
	    dDa = (c1 - c1 * c2 * c2 - 2 * c2) * dnormc1 * imahalSquare * imahal *
	    ifrech1Square * ifrech2 + (c2 - c1 * c1 *c2 - 2 * c1) * dnormc2 *
	    imahalSquare * imahal * ifrech1 * ifrech2Square,
	    jacCommonSigma = weights[currentPair] * (dAa + (dBa * C + B * dCa + dDa) / (B*C + D));

	  hess[k * nPairs + currentPair] = -(*cov12 * distVec[nPairs + currentPair] - *cov22 * distVec[currentPair]) *
	    (*cov12 * distVec[nPairs + currentPair] - *cov22 * distVec[currentPair]) /
	    (2 * det * det * mahalDist[currentPair]) * jacCommonSigma;

	  hess[(*nObs + k) * nPairs + currentPair] = (*cov11 * distVec[nPairs + currentPair] -
						      *cov12 * distVec[currentPair]) *
	    (*cov12 * distVec[nPairs + currentPair] - *cov22 * distVec[currentPair]) /
	    (det * det * mahalDist[currentPair]) * jacCommonSigma;
	  hess[(2 * *nObs + k) * nPairs + currentPair] = -(*cov11 * distVec[nPairs + currentPair] -
							   *cov12 * distVec[currentPair]) *
	    (*cov11 * distVec[nPairs + currentPair] -  *cov12 * distVec[currentPair]) /
	    (2 * det * det * mahalDist[currentPair]) * jacCommonSigma;

	  grad[k] += hess[k * nPairs + currentPair];
	  grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];
	  grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];
	}
      }
    }
  }

  if (*fitmarge){
    int start = 3;
    marginalPartSmith(&start, nObs, nSite, data, frech, mahalDist, locs, scales, shapes,
		      trendlocs, trendscales, trendshapes, nloccoeff, nscalecoeff, nshapecoeff,
		      ntemploccoeff, ntempscalecoeff, ntempshapecoeff, locdsgnmat, scaledsgnmat,
		      shapedsgnmat, tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, weights,
		      hess, grad);
  }

  return;
}

void smithstderr3d(double *data, double *distVec, int *nSite, int *nObs, double *locdsgnmat,
		   int *nloccoeff, double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		   int *nshapecoeff, double *tempdsgnmatloc, int *ntemploccoeff,
		   double *tempdsgnmatscale, int *ntempscalecoeff, double *tempdsgnmatshape,
		   int *ntempshapecoeff, double *loccoeff, double *scalecoeff, double *shapecoeff,
		   double *temploccoeff, double *tempscalecoeff, double *tempshapecoeff, double *cov11,
		   double *cov12, double *cov13, double *cov22, double *cov23, double *cov33, int *fitmarge,
		   int *usetempcov, double *weights, double *hess, double *grad){

  //This is the Smith model. It computes the hessian of the pairwise log-likelihood

  const int nPairs = *nSite * (*nSite - 1) / 2,
    flag = usetempcov[0] + usetempcov[1] + usetempcov[2],
    det = *cov11 * *cov22 * *cov33 - *cov12 * *cov12 * *cov33 -
    *cov11 * *cov23 * *cov23 + 2 * *cov12 * *cov13 * *cov23 -
    *cov13 * *cov13 * *cov22;

  int i, currentPair = -1;
  double *mahalDist, *locs, *scales, *shapes, *trendlocs, *trendscales, *trendshapes, *jac,
    *frech;

  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  trendlocs = (double *)R_alloc(*nObs, sizeof(double));
  trendscales = (double *)R_alloc(*nObs, sizeof(double));
  trendshapes = (double *)R_alloc(*nObs, sizeof(double));

  /* We have to intialized them to 0 as we don't know if
     dsgnmat2temptrend will be called */
  for (i=*nObs;i--;)
    trendlocs[i] = trendscales[i] = trendshapes[i] = 0;

  //Computing the Mahalanobis distance
  mahalDistFct3d(distVec, nPairs, cov11, cov12, cov13, cov22, cov23, cov33, mahalDist);

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){

    dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff, shapecoeff, *nSite,
		  *nloccoeff, *nscalecoeff, *nshapecoeff, locs, scales, shapes);

    //Stage 1: Transformation to unit Frechet
    if (flag) {
      dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, temploccoeff,
			tempscalecoeff, tempshapecoeff, *nSite, *nObs, usetempcov, *ntemploccoeff,
			*ntempscalecoeff, *ntempshapecoeff, trendlocs, trendscales, trendshapes);

      gev2frechTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs, trendscales,
		     trendshapes, jac, frech);
    }

    else
      gev2frech(data, *nObs, *nSite, locs, scales, shapes, jac, frech);
  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];

  //Stage 2: Hessian computations
  // a- Covariance matrix part
  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){
      currentPair++;

      int k;
      double imahal = 1 / mahalDist[currentPair], imahalSquare = imahal * imahal;

      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    hess[k * nPairs + currentPair] = hess[(*nObs + k) * nPairs + currentPair] =
	      hess[(2 * *nObs + k) * nPairs + currentPair] = hess[(3 * *nObs + k) * nPairs + currentPair] =
	      hess[(4 * *nObs + k) * nPairs + currentPair] = hess[(5 * *nObs + k) * nPairs + currentPair] =
	      NA_REAL;

	    continue;
	  }

	  double ifrech1 = 1 / frech[k + i * *nObs], ifrech2 = 1 / frech[k + j * *nObs],
	    ifrech1Square = ifrech1 * ifrech1, ifrech2Square = ifrech2 * ifrech2,
	    c1 = log(frech[k + j * *nObs] * ifrech1) * imahal + 0.5 * mahalDist[currentPair],
	    c2 = mahalDist[currentPair] - c1,
	    dnormc1 = dnorm(c1, 0., 1., 0), pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0), pnormc2 = pnorm(c2, 0., 1., 1, 0),
	    //A = - pnormc1 * ifrech1 - pnormc2 * ifrech2,
	    B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	    dnormc2 * imahal * ifrech2Square,
	    C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	    dnormc1 * imahal * ifrech1Square,
	    D = c2 * dnormc1 * imahalSquare * ifrech1Square * ifrech2 +
	    c1 * dnormc2 * imahalSquare * ifrech1 * ifrech2Square,
	    dAa = - c2 * dnormc1 * imahal * ifrech1 - c1 * dnormc2 * imahal * ifrech2,
	    dBa = (c1 * c1 - 1) * dnormc2 * imahalSquare * ifrech2Square +
	    (1 + c1 * c2 ) * dnormc1 * imahalSquare * ifrech1 * ifrech2,
	    dCa = (c2 * c2 - 1) * dnormc1 * imahalSquare * ifrech1Square +
	    (1 + c1 * c2) * dnormc2 * imahalSquare * ifrech1 * ifrech2,
	    dDa = (c1 - c1 * c2 * c2 - 2 * c2) * dnormc1 * imahal * imahalSquare *
	    ifrech1Square * ifrech2 + (c2 - c1 * c1 *c2 - 2 * c1) * dnormc2 *
	    imahalSquare * imahal * ifrech1 * ifrech2Square,
	    jacCommonSigma = weights[currentPair] * (dAa + (dBa * C + B * dCa + dDa) / (B*C + D));

	  hess[k * nPairs + currentPair] = -((*cov22 * *cov33 - *cov23 * *cov23) * distVec[currentPair] +
					     (*cov13 * *cov23 - *cov12 * *cov33) * distVec[nPairs + currentPair] +
					     (*cov12 * *cov23 - *cov13 * *cov22) * distVec[2 * nPairs + currentPair]) *
	    ((*cov22 * *cov33 - *cov23 * *cov23) * distVec[currentPair] +
	     (*cov13 * *cov23 - *cov12 * *cov33) * distVec[nPairs + currentPair] +
	     (*cov12 * *cov23 - *cov13 * *cov22) * distVec[2 * nPairs + currentPair]) /
	    (2 * det * det * mahalDist[currentPair]) * jacCommonSigma;

	  hess[(*nObs + k) * nPairs + currentPair] = ((*cov12 * *cov33 - *cov13 * *cov23) * distVec[currentPair] +
						      (*cov13 * *cov13 - *cov11 * *cov33) * distVec[nPairs + currentPair] +
						      (*cov11 * *cov23 - *cov12 * *cov13) * distVec[2 * nPairs + currentPair]) *
	    ((*cov12 * *cov23 - *cov13 * *cov22) * distVec[2 * nPairs + currentPair] +
	     (*cov22 * *cov33 - *cov23 * *cov23) * distVec[currentPair] +
	     (*cov13 * *cov23 - *cov12 * *cov33) * distVec[nPairs + currentPair]) /
	    (det * det * mahalDist[currentPair]) * jacCommonSigma;

	  hess[(2 * *nObs + k) * nPairs + currentPair] = -((*cov22 * *cov33 - *cov23 * *cov23) * distVec[currentPair] +
							   (*cov13 * *cov23 - *cov12 * *cov33) * distVec[nPairs + currentPair] +
							   (*cov12 * *cov23 - *cov13 * *cov22) * distVec[2 * nPairs + currentPair]) *
	    ((*cov11 * *cov22 - *cov12 * *cov12) * distVec[2 * nPairs + currentPair] +
	     (*cov12 * *cov23 - *cov13 * *cov22) * distVec[currentPair] +
	     (*cov12 * *cov13 - *cov11 * *cov23) * distVec[nPairs + currentPair]) /
	    (det * det * mahalDist[currentPair]) * jacCommonSigma;

	  hess[(3 * *nObs + k) * nPairs + currentPair] = -
	    ((*cov11 * *cov23 - *cov12 * *cov13) * distVec[2 * nPairs + currentPair] +
	     (*cov13 * *cov13 - *cov11 * *cov33) * distVec[nPairs + currentPair] +
	     (*cov12 * *cov33 - *cov13 * *cov23) * distVec[currentPair]) *
	    ((*cov11 * *cov23 - *cov12 * *cov13) * distVec[2 * nPairs + currentPair] +
	     (*cov13 * *cov13 - *cov11 * *cov33) * distVec[nPairs + currentPair] +
	     (*cov12 * *cov33 - *cov13 * *cov23) * distVec[currentPair]) /
	    (2 * det * det * mahalDist[currentPair]) * jacCommonSigma;

	  hess[(4 * *nObs + k) * nPairs + currentPair] =
	    ((*cov11 * *cov23 - *cov12 * *cov13) * distVec[2 * nPairs + currentPair] +
	     (*cov12 * *cov33 - *cov13 * *cov23) * distVec[currentPair] +
	     (*cov13 * *cov13 - *cov11 * *cov33) * distVec[nPairs + currentPair]) *
	    ((*cov11 * *cov22 - *cov12 * *cov12) * distVec[2 * nPairs + currentPair] +
	     (*cov12 * *cov13 - *cov11 * *cov23) * distVec[nPairs + currentPair] +
	     (*cov12 * *cov23 - *cov13 * *cov22) * distVec[currentPair]) /
	    (det * det * mahalDist[currentPair]) * jacCommonSigma;

	  hess[(5 * *nObs + k) * nPairs + currentPair] = -
	    ((*cov11 * *cov22 - *cov12 * *cov12) * distVec[2 * nPairs + currentPair] +
	     (*cov12 * *cov13 - *cov11 * *cov23) * distVec[nPairs + currentPair] +
	     (*cov12 * *cov23 - *cov13 * *cov22) * distVec[currentPair]) *
	    ((*cov11 * *cov22 - *cov12 * *cov12) * distVec[2 * nPairs + currentPair] +
	     (*cov12 * *cov13 - *cov11 * *cov23) * distVec[nPairs + currentPair] +
	     (*cov12 * *cov23 - *cov13 * *cov22) * distVec[currentPair]) /
	    (2 * det * det * mahalDist[currentPair]) * jacCommonSigma;

	  if (!ISNA(hess[k * nPairs + currentPair])){
	    grad[k] += hess[k * nPairs + currentPair];
	    grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];
	    grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];
	    grad[3 * *nObs + k] += hess[(3 * *nObs + k) * nPairs + currentPair];
	    grad[4 * *nObs + k] += hess[(4 * *nObs + k) * nPairs + currentPair];
	    grad[5 * *nObs + k] += hess[(5 * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }

  if (*fitmarge){
    int start = 6;
    marginalPartSmith(&start, nObs, nSite, data, frech, mahalDist, locs, scales, shapes,
		      trendlocs, trendscales, trendshapes, nloccoeff, nscalecoeff, nshapecoeff,
		      ntemploccoeff, ntempscalecoeff, ntempshapecoeff, locdsgnmat, scaledsgnmat,
		      shapedsgnmat, tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, weights, hess, grad);
  }

  return;
}


void schlatherstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		     double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		     double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
		     int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
		     double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
		     double *scalecoeff, double *shapecoeff, double *temploccoeff,
		     double *tempscalecoeff, double *tempshapecoeff, double *nugget, double *range,
		     double *smooth, double *smooth2, int *fitmarge, int *usetempcov, double *weights,
		     double *hess, double *grad){

  /* This is the Schlather model. It computes the hessian of the
     pairwise log-likelihood */

  const int nPairs = *nSite * (*nSite - 1) / 2,
    flag = usetempcov[0] + usetempcov[1] + usetempcov[2];
  int i, currentPair = -1, nCorPar = 3;
  double *rho, *locs, *scales, *shapes, *trendlocs, *trendscales, *trendshapes, *jac, *frech, h;
  // h is used for numerical derivative (whitmat and bessel)
  const double sill = 1 - *nugget;

  h = sqrt(EPS) * *smooth;
  double temp = *smooth + h;
  h = temp - *smooth;

  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  trendlocs = (double *)R_alloc(*nObs, sizeof(double));
  trendscales = (double *)R_alloc(*nObs, sizeof(double));
  trendshapes = (double *)R_alloc(*nObs, sizeof(double));

  /* We have to intialized them to 0 as we don't know if
     dsgnmat2temptrend will be called */
  for (i=*nObs;i--;)
    trendlocs[i] = trendscales[i] = trendshapes[i] = 0;

  //Stage 0: Compute the covariance at each location
  switch (*covmod){
  case 1:
    whittleMatern(dist, nPairs, *nugget, sill, *range, *smooth, rho);
    break;
  case 2:
    cauchy(dist, nPairs, *nugget, sill, *range, *smooth, rho);
    break;
  case 3:
    powerExp(dist, nPairs, *nugget, sill, *range, *smooth, rho);
    break;
  case 4:
    //Here we use 0 for dim as we don't care for the computation of
    //the hessian
    bessel(dist, nPairs, 0, *nugget, sill, *range, *smooth, rho);
    break;
  case 5:
    caugen(dist, nPairs, *nugget, sill, *range, *smooth, *smooth2, rho);
  }

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){

    dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff,
		  shapecoeff, *nSite, *nloccoeff, *nscalecoeff, *nshapecoeff,
		  locs, scales, shapes);

    //Stage 1: Transformation to unit Frechet
    if (flag) {
      dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, temploccoeff,
			tempscalecoeff, tempshapecoeff, *nSite, *nObs, usetempcov, *ntemploccoeff,
			*ntempscalecoeff, *ntempshapecoeff, trendlocs, trendscales, trendshapes);

      gev2frechTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs, trendscales,
		     trendshapes, jac, frech);
    }

    else
      gev2frech(data, *nObs, *nSite, locs, scales, shapes, jac, frech);

  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];

  //Stage 2: Hessian computations;
  // a- Covariance part
  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){

      currentPair++;

      int k;
      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    hess[k * nPairs + currentPair] = hess[(*nObs + k) * nPairs + currentPair] =
	      hess[(2 * *nObs + k) * nPairs + currentPair] = NA_REAL;

	    if (*covmod == 5)
	      hess[(3 * *nObs + k) * nPairs + currentPair] = NA_REAL;

	    continue;
	  }

	  double c1 = sqrt(frech[k + i * *nObs] * frech[k + i * *nObs] + frech[k + j * *nObs] * frech[k + j * *nObs] -
			   2 * frech[k + i * *nObs] * frech[k + j * *nObs] * rho[currentPair]),
	    B = (1 - rho[currentPair] * rho[currentPair]) / (2 * c1 * c1 * c1),
	    C = (- rho[currentPair] * frech[k + i * *nObs] + c1 + frech[k + j * *nObs]) /
	    (2 * c1 * frech[k + i * *nObs] * frech[k + i * *nObs]),
	    D = (-rho[currentPair] * frech[k + j * *nObs] + c1 + frech[k + i * *nObs]) /
	    (2 * c1 * frech[k + j * *nObs] * frech[k + j * *nObs]),
	    dArho =  1 / (2 * c1),
	    dBrho = - rho[currentPair] / (c1 * c1 * c1) + 3 * (1 - rho[currentPair] * rho[currentPair]) *
	    frech[k + i * *nObs] * frech[k + j * *nObs] / (2 * c1 * c1 * c1 * c1 * c1),
	    dCrho = (-frech[k + i * *nObs] + frech[k + j * *nObs] * rho[currentPair]) / (2 * c1 * c1 * c1),
	    dDrho = (-frech[k + j * *nObs] + frech[k + i * *nObs] * rho[currentPair]) / (2 * c1 * c1 * c1),
	    jacCommonRho = weights[currentPair] * (dArho + (dBrho + dCrho * D + C * dDrho) / (B + C * D));

	  hess[k * nPairs + currentPair] = -rho[currentPair] / sill * jacCommonRho;
	  grad[k] += hess[k * nPairs + currentPair];

	  switch (*covmod){
	  case 1:
	    //i.e. Whittle-Matern
	    hess[(*nObs + k) * nPairs + currentPair] = rho[currentPair] *
	      (-2 * *smooth / *range + dist[currentPair] *
	       bessel_k(dist[currentPair] / *range, *smooth + 1, 1) /
	       (bessel_k(dist[currentPair] / *range, *smooth, 1) * *range * *range)) *
	      jacCommonRho;
	    /* There's no closed form for the partial derivative of the
	       BesselK function w.r.t. smooth. We use finite differences
	       for this... */
	    hess[(2 * *nObs + k) * nPairs + currentPair] = rho[currentPair] *
	      (-M_LN2 - digamma(*smooth) + log(dist[currentPair] / *range) +
	       (bessel_k(dist[currentPair] / *range, *smooth + h, 1) -
		bessel_k(dist[currentPair] / *range, *smooth - h, 1)) /
	       (2 * h * bessel_k(dist[currentPair] / *range, *smooth, 1))) *
	      jacCommonRho;
	    break;
	  case 2:
	    //i.e. cauchy
	    hess[(*nObs + k) * nPairs + currentPair] = 2 * dist[currentPair] * dist[currentPair] *
	      sill * *smooth / (*range * *range * *range) *
	      R_pow(1 + dist[currentPair] * dist[currentPair] / (*range * *range),
		    - *smooth - 1) * jacCommonRho;
	    hess[(2 * *nObs + k) * nPairs + currentPair] = -rho[currentPair] *
	      log1p(dist[currentPair] *dist[currentPair] / ( *range * *range)) *
	      jacCommonRho;
	    break;
	  case 3:
	    //i.e. powered exponential
	    hess[(*nObs + k) * nPairs + currentPair] = rho[currentPair] * *smooth / *range *
	      R_pow(dist[currentPair] / *range, *smooth) * jacCommonRho;
	    hess[(2 * *nObs + k) * nPairs + currentPair] = -rho[currentPair] *
	      R_pow(dist[currentPair] / *range, *smooth) * log(dist[currentPair] / *range) *
	      jacCommonRho;
	    break;
	  case 4:
	    //i.e. Bessel
	    hess[(*nObs + k) * nPairs + currentPair] = sill *
	      R_pow(2 * *range / dist[currentPair], *smooth) *
	      dist[currentPair] / (*range * *range) * gammafn(*smooth + 1) *
	      bessel_j(dist[currentPair] / *range, *smooth + 1) * jacCommonRho;
	    /* There's no closed form for the partial derivative of the
	       BesselJ function w.r.t. smooth. We use finite differences
	       for this... */
	    hess[(2 * *nObs + k) * nPairs + currentPair] = rho[currentPair] *
	      (M_LN2 - log(dist[currentPair] / *range) + digamma(*smooth + 1) +
	       (bessel_j(dist[currentPair] / *range, *smooth + h) -
		bessel_j(dist[currentPair] / *range, *smooth)) /
	       (h * bessel_j(dist[currentPair] / *range, *smooth))) *
	      jacCommonRho;
	    break;
	  case 5:
	    //i.e. Generalized Cauchy
	    hess[(*nObs + k) * nPairs + currentPair] = rho[currentPair] * *smooth / *range /
	      (1 + R_pow(*range / dist[currentPair], *smooth2)) *
	      jacCommonRho;
	    hess[(2 * *nObs + k) * nPairs + currentPair] = -log1p(R_pow(dist[currentPair] / *range, *smooth2)) *
	      rho[currentPair] / *smooth2 *jacCommonRho;
	    hess[(3 * *nObs + k) * nPairs + currentPair] = (log1p(R_pow(dist[currentPair] / *range, *smooth2)) /
							    *smooth2 - log(dist[currentPair] / *range) /
							    (1 + R_pow(*range / dist[currentPair], *smooth2))) *
	      *smooth / *smooth2 * rho[currentPair] * jacCommonRho;
	    grad[3 * *nObs + k] += hess[(3 * *nObs + k) * nPairs + currentPair];
	    //The caugen has 2 smooth paramaters so add + 1 to nCorPar.
	    nCorPar = 4;
	    break;
	  }

	  grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];
	  grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];
	}
      }
    }
  }

  if (*fitmarge)
    marginalPartSchlat(&nCorPar, nObs, nSite, data, frech, rho, locs, scales, shapes, trendlocs,
		       trendscales, trendshapes, nloccoeff, nscalecoeff, nshapecoeff,
		       ntemploccoeff, ntempscalecoeff, ntempshapecoeff, locdsgnmat, scaledsgnmat,
		       shapedsgnmat, tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, weights, hess,
		       grad);

  return;
}

void schlatherindstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
			double *locdsgnmat, int *nloccoeff, double *scaledsgnmat,
			int *nscalecoeff, double *shapedsgnmat,	int *nshapecoeff,
			double *tempdsgnmatloc, int *ntemploccoeff, double *tempdsgnmatscale,
			int *ntempscalecoeff, double *tempdsgnmatshape, int *ntempshapecoeff,
			double *loccoeff, double *scalecoeff, double *shapecoeff,
			double *temploccoeff, double *tempscalecoeff, double *tempshapecoeff,
			double *alpha, double *nugget, double *range, double *smooth,
			double *smooth2, int *fitmarge, int *usetempcov, double *weights, double *hess,
			double *grad){

  /* This is the independent Schlather model. It computes the hessian
     of the pairwise log-likelihood */

  const int nPairs = *nSite * (*nSite - 1) / 2,
    flag = usetempcov[0] + usetempcov[1] + usetempcov[2];
  int i, currentPair = -1, nCorPar = 4;
  double *rho, *locs, *scales, *shapes, *trendlocs, *trendscales, *trendshapes, *jac, *frech, h;
  const double sill = 1 - *nugget;

  h = sqrt(EPS) * *smooth;
  double temp = *smooth + h;
  h = temp - *smooth;

  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  trendlocs = (double *)R_alloc(*nObs, sizeof(double));
  trendscales = (double *)R_alloc(*nObs, sizeof(double));
  trendshapes = (double *)R_alloc(*nObs, sizeof(double));

  /* We have to intialized them to 0 as we don't know if
     dsgnmat2temptrend will be called */
  for (i=*nObs;i--;)
    trendlocs[i] = trendscales[i] = trendshapes[i] = 0;

  //Stage 0: Compute the covariance at each location
  switch (*covmod){
  case 1:
    whittleMatern(dist, nPairs, *nugget, sill, *range, *smooth, rho);
    break;
  case 2:
    cauchy(dist, nPairs, *nugget, sill, *range, *smooth, rho);
    break;
  case 3:
    powerExp(dist, nPairs, *nugget, sill, *range, *smooth, rho);
    break;
  case 4:
    //Here we use 0 for dim as we don't care for the computation of
    //the hessian
    bessel(dist, nPairs, 0, *nugget, sill, *range, *smooth, rho);
    break;
  case 5:
    caugen(dist, nPairs, *nugget, sill, *range, *smooth, *smooth2, rho);
    break;
  }

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){

    dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff, shapecoeff,
		  *nSite, *nloccoeff, *nscalecoeff, *nshapecoeff, locs, scales, shapes);

    //Stage 1: Transformation to unit Frechet
    if (flag) {
      dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, temploccoeff,
			tempscalecoeff, tempshapecoeff, *nSite, *nObs, usetempcov, *ntemploccoeff,
			*ntempscalecoeff, *ntempshapecoeff, trendlocs, trendscales, trendshapes);

      gev2frechTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs, trendscales,
		     trendshapes, jac, frech);
    }

    else
      gev2frech(data, *nObs, *nSite, locs, scales, shapes, jac, frech);

  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];

  //Stage 2: Hessian computations;
  // a- Covariance part
  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){

      currentPair++;
      int k;

      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    hess[k * nPairs + currentPair] = hess[(*nObs + k) * nPairs + currentPair] =
	      hess[(2 * *nObs + k) * nPairs + currentPair] = hess[(3 * *nObs + k) * nPairs + currentPair] =
	      NA_REAL;

	    if (*covmod == 5)
	      hess[(4 * *nObs + k) * nPairs + currentPair] = NA_REAL;

	    continue;
	  }

	  double frech1Square = frech[k + i * *nObs] * frech[k + i * *nObs],
	    frech2Square = frech[k + j * *nObs] * frech[k + j * *nObs],
	    c1 = sqrt(frech1Square + frech2Square - 2 * frech[k + i * *nObs] *
		      frech[k + j * *nObs] * rho[currentPair]),
	    B = (1 - *alpha) * (1 - rho[currentPair] * rho[currentPair]) / (2 * c1 * c1 * c1),
	    C = (*alpha - 1) * (rho[currentPair] * frech[k + i * *nObs] - c1 -
				frech[k + j * *nObs]) / (2 * c1 * frech[k + i * *nObs] *
							 frech[k + i * *nObs]) +
	    *alpha / (frech[k + i * *nObs] * frech[k + i * *nObs]),
	    D = (*alpha - 1) * (rho[currentPair] * frech[k + j * *nObs] - c1 -
				frech[k + i * *nObs]) / (2 * c1 * frech[k + j * *nObs] *
							 frech[k + j * *nObs]) +
	    *alpha / (frech[k + j * *nObs] * frech[k + j * *nObs]),
	    dAalpha = (c1 - frech[k + i * *nObs] - frech[k + j * *nObs]) /
	    (2 * frech[k + i * *nObs] * frech[k + j * *nObs]),
	    dArho =  (1 - *alpha) / (2 * c1),
	    dBalpha = (rho[currentPair] * rho[currentPair] - 1) / (2 * c1 * c1 * c1),
	    dBrho = (3 * (1 - rho[currentPair] * rho[currentPair]) * frech[k + i * *nObs] *
		     frech[k + j * *nObs] / (2 * c1 * c1 * c1 * c1 * c1) - rho[currentPair] /
		     (c1 * c1 * c1)) * (1 - *alpha),
	    dCalpha = (rho[currentPair] * frech[k + i * *nObs] + c1 - frech[k + j * *nObs]) /
	    (2 * c1 * frech[k + i * *nObs] * frech[k + i * *nObs]),
	    dCrho = (*alpha - 1) * (frech[k + i * *nObs] - frech[k + j * *nObs] *
				    rho[currentPair]) / (2 * c1 * c1 * c1),
	    dDalpha = (rho[currentPair] * frech[k + j * *nObs] + c1 - frech[k + i * *nObs]) /
	    (2 * c1 * frech[k + j * *nObs] * frech[k + j * *nObs]),
	    dDrho = (*alpha - 1) * (frech[k + j * *nObs] - frech[k + i * *nObs] *
				    rho[currentPair]) / (2 * c1 * c1 * c1),
	    jacCommonRho = weights[currentPair] * (dArho + (dBrho + dCrho * D + dDrho * C) / (B + C * D));

	  hess[k * nPairs + currentPair] = weights[currentPair] *
	    (dAalpha + (dBalpha + dCalpha * D + dDalpha * C) / (B + C * D));
	  hess[(*nObs + k) * nPairs + currentPair] = -rho[currentPair] / sill * jacCommonRho;

	  grad[k] += hess[k * nPairs + currentPair];
	  grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];

	  switch (*covmod){
	  case 1:
	    //i.e. Whittle-Matern
	    hess[(2 * *nObs + k) * nPairs + currentPair] = rho[currentPair] *
	      (-2 * *smooth / *range + dist[currentPair] *
	       bessel_k(dist[currentPair] / *range, *smooth + 1, 1) /
	       bessel_k(dist[currentPair] / *range, *smooth, 1) * *range * *range) *
	      jacCommonRho;
	    /* There's no closed form for the partial derivative of the
	       BesselK function w.r.t. smooth. We use finite differences
	       for this... */
	    hess[(3 * *nObs + k) * nPairs + currentPair] = rho[currentPair] *
	      (-M_LN2 - digamma(*smooth) + log(dist[currentPair] / *range) +
	       (bessel_k(dist[currentPair] / *range, *smooth + h, 1) -
		bessel_k(dist[currentPair] / *range, *smooth - h, 1)) /
	       (2 * h * bessel_k(dist[currentPair] / *range, *smooth, 1))) *
	      jacCommonRho;
	    break;
	  case 2:
	    //i.e. cauchy
	    hess[(2 * *nObs + k) * nPairs + currentPair] = 2 * dist[currentPair] * dist[currentPair] *
	      sill * *smooth / (*range * *range * *range) *
	      R_pow(1 + dist[currentPair] * dist[currentPair] / (*range * *range),
		    - *smooth - 1) * jacCommonRho;
	    hess[(3 * *nObs + k) * nPairs + currentPair] = -rho[currentPair] *
	      log1p(dist[currentPair] * dist[currentPair] / (*range * *range)) *
	      jacCommonRho;
	    break;
	  case 3:
	    //i.e. powered exponential
	    hess[(2 * *nObs + k) * nPairs + currentPair] = rho[currentPair] * *smooth /
	      *range * R_pow(dist[currentPair] / *range, *smooth) * jacCommonRho;
	    hess[(3 * *nObs + k) * nPairs + currentPair] = -rho[currentPair] *
	      R_pow(dist[currentPair] / *range, *smooth) * log(dist[currentPair] / *range) *
	      jacCommonRho;
	    break;
	  case 4:
	    //i.e. Bessel
	    hess[(2 * *nObs + k) * nPairs + currentPair] = sill *
	      R_pow(2 * *range / dist[currentPair], *smooth) *
	      dist[currentPair] / (*range * *range) * gammafn(*smooth + 1) *
	      bessel_j(dist[currentPair] / *range, *smooth + 1) * jacCommonRho;
	    /* There's no closed form for the partial derivative of the
	       BesselJ function w.r.t. smooth. We use finite differences
	       for this... */
	    hess[(3 * *nObs + k) * nPairs + currentPair] = rho[currentPair] *
	      (M_LN2 - log(dist[currentPair] / *range) + digamma(*smooth + 1) +
	       (bessel_j(dist[currentPair] / *range, *smooth + h) -
		bessel_j(dist[currentPair] / *range, *smooth- h)) /
	       (2 * h *  - bessel_j(dist[currentPair] / *range, *smooth))) *
	      jacCommonRho;
	    break;
	  case 5:
	    //i.e. Generalized Cauchy
	    hess[(2 * *nObs + k) * nPairs + currentPair] = rho[currentPair] * *smooth / *range /
	      (1 + R_pow(*range / dist[currentPair], *smooth2)) *
	      jacCommonRho;
	    hess[(3 * *nObs + k) * nPairs + currentPair] = -log1p(R_pow(dist[currentPair] / *range, *smooth2)) *
	      rho[currentPair] / *smooth2 *jacCommonRho;
	    hess[(4 * *nObs + k) * nPairs + currentPair] = (log1p(R_pow(dist[currentPair] / *range, *smooth2)) /
							    *smooth2 - log(dist[currentPair] / *range) /
							    (1 + R_pow(*range / dist[currentPair], *smooth2))) *
	      *smooth / *smooth2 * rho[currentPair] * jacCommonRho;
	    grad[4 * *nObs + k] += hess[(4 * *nObs + k) * nPairs + currentPair];
	    //The caugen has 2 smooth paramaters so add + 1 to nCorPar.
	    nCorPar = 5;
	    break;
	  }
	  grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];
	  grad[3 * *nObs + k] += hess[(3 * *nObs + k) * nPairs + currentPair];
	}
      }
    }
  }

  if (*fitmarge)
    marginalPartiSchlat(&nCorPar, nObs, nSite, data, frech, alpha, rho, locs, scales, shapes,
			trendlocs, trendscales, trendshapes, nloccoeff, nscalecoeff, nshapecoeff,
			ntemploccoeff, ntempscalecoeff, ntempshapecoeff, locdsgnmat, scaledsgnmat,
			shapedsgnmat, tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, weights, hess,
			grad);

  return;
}

void geomgaussstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		     double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		     double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
		     int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
		     double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
		     double *scalecoeff, double *shapecoeff, double *temploccoeff,
		     double *tempscalecoeff, double *tempshapecoeff, double *sigma2,
		     double *nugget, double *range, double *smooth, double *smooth2,
		     int *fitmarge, int *usetempcov, double *weights, double *hess, double *grad){
  /* This function computes the hessian of the log-pairwise
     likelihood for the geometric gaussian model.

     Remember that this model shares the same bivariate density as the
     Smith model except that the Mahalanobis distance is modified. */

  const int nPairs = *nSite * (*nSite - 1) / 2,
    flag = usetempcov[0] + usetempcov[1] + usetempcov[2];
  int i, currentPair = -1, nCorPar = 4;
  double *mahalDist, *locs, *scales, *shapes, *jac, *frech, *trendlocs,
    *trendscales, *trendshapes, h;
  const double sill = 1 - *nugget;

  h = sqrt(EPS) * *smooth;
  double temp = *smooth + h;
  h = temp - *smooth;

  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  trendlocs = (double *)R_alloc(*nObs, sizeof(double));
  trendscales = (double *)R_alloc(*nObs, sizeof(double));
  trendshapes = (double *)R_alloc(*nObs, sizeof(double));

  /* We have to intialized them to 0 as we don't know if
     dsgnmat2temptrend will be called */
  for (i=*nObs;i--;)
    trendlocs[i] = trendscales[i] = trendshapes[i] = 0;

  //Stage 0: Compute the covariance at each location
  geomCovariance(dist, nPairs, 0, *covmod, *sigma2, 1.0e6, *nugget, *range, *smooth, *smooth2,
		 mahalDist);

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){

    dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff, shapecoeff,
		  *nSite, *nloccoeff, *nscalecoeff, *nshapecoeff, locs, scales, shapes);

    //Stage 1: Transformation to unit Frechet
    if (flag) {
      dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, temploccoeff,
			tempscalecoeff, tempshapecoeff, *nSite, *nObs, usetempcov, *ntemploccoeff,
			*ntempscalecoeff, *ntempshapecoeff, trendlocs, trendscales, trendshapes);

      gev2frechTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs, trendscales,
		     trendshapes, jac, frech);
    }

    else
      gev2frech(data, *nObs, *nSite, locs, scales, shapes, jac, frech);

  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];

  //Stage 2: Hessian computations;
  // a- Covariance part
  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){
      currentPair++;

      int k;
      double imahal = 1 / mahalDist[currentPair], imahalSquare = imahal * imahal,
	rho = 1 - mahalDist[currentPair] * mahalDist[currentPair] / (2 * *sigma2);

      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    hess[k * nPairs + currentPair] = hess[(*nObs + k) * nPairs + currentPair] =
	      hess[(2 * *nObs + k) * nPairs + currentPair] = hess[(3 * *nObs + k) * nPairs + currentPair] =
	      NA_REAL;

	    if (*covmod == 5)
	      hess[(4 * *nObs + k) * nPairs + currentPair] = NA_REAL;

	    continue;
	  }

	  double ifrech1 = 1 / frech[k + i * *nObs], ifrech2 = 1 / frech[k + j * *nObs],
	    ifrech1Square = ifrech1 * ifrech1, ifrech2Square = ifrech2 * ifrech2,
	    c1 = log(frech[k + j * *nObs] * ifrech1) * imahal + 0.5 * mahalDist[currentPair],
	    c2 = mahalDist[currentPair] - c1,
	    dnormc1 = dnorm(c1, 0., 1., 0), pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0), pnormc2 = pnorm(c2, 0., 1., 1, 0),
	    //A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	    B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	    dnormc2 * imahal * ifrech2Square,
	    C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	    dnormc1 * imahal * ifrech1Square,
	    D = c2 * dnormc1 * imahalSquare * ifrech1Square * ifrech2 +
	    c1 * dnormc2 * imahalSquare * ifrech1 * ifrech2Square,
	    dAa = - c2 * dnormc1 * imahal * ifrech1 - c1 * dnormc2 * imahal * ifrech2,
	    dBa = (c1 * c1 - 1) * dnormc2 * imahalSquare * ifrech2Square +
	    (1 + c1 * c2 ) * dnormc1 * imahalSquare * ifrech1 * ifrech2,
	    dCa = (c2 * c2 - 1) * dnormc1 * imahalSquare * ifrech1Square +
	    (1 + c1 * c2) * dnormc2 * imahalSquare * ifrech1 * ifrech2,
	    dDa = (c1 - c1 * c2 * c2 - 2 * c2) * dnormc1 * imahalSquare * imahal *
	    ifrech1Square * ifrech2 + (c2 - c1 * c1 *c2 - 2 * c1) * dnormc2 *
	    imahalSquare * imahal * ifrech1 * ifrech2Square,
	    jacCommon = weights[currentPair] * (dAa + (dBa * C + B * dCa + dDa) / (B*C + D));

	  hess[k * nPairs + currentPair] = mahalDist[currentPair] / (2 * *sigma2) * jacCommon;
	  hess[(*nObs + k) * nPairs + currentPair] = *sigma2 * rho / (mahalDist[currentPair] * sill) *
	    jacCommon;

	  grad[k] += hess[k * nPairs + currentPair];
	  grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];

	  switch (*covmod){
	  case 1:
	    //i.e. Whittle-Matern
	    hess[(2 * *nObs + k) * nPairs + currentPair] = *sigma2 * rho / mahalDist[currentPair] *
	      (2 * *smooth / *range - dist[currentPair] *
	       bessel_k(dist[currentPair] / *range, *smooth + 1, 1) /
	       (bessel_k(dist[currentPair] / *range, *smooth, 1) * *range * *range)) *
	      jacCommon;
	    /* There's no closed form for the partial derivative of the
	       BesselK function w.r.t. smooth. We use finite differences
	       for this... */
	    hess[(3 * *nObs + k) * nPairs + currentPair] = -*sigma2 * rho / mahalDist[currentPair] *
	      (-M_LN2 - digamma(*smooth) + log(dist[currentPair] / *range) +
	       (bessel_k(dist[currentPair] / *range, *smooth + h, 1) -
		bessel_k(dist[currentPair] / *range, *smooth - h, 1)) /
	       (2 * h * bessel_k(dist[currentPair] / *range, *smooth, 1))) *
	      jacCommon;
	    break;
	  case 2:
	    //i.e. cauchy
	    hess[(2 * *nObs + k) * nPairs + currentPair] = -2 * *sigma2 * dist[currentPair] * dist[currentPair] *
	      sill * *smooth * R_pow(1 + dist[currentPair] * dist[currentPair] /
				     (*range * *range), - *smooth - 1) /
	      (*range * *range * *range * mahalDist[currentPair]) * jacCommon;
	    hess[(3 * *nObs + k) * nPairs + currentPair] = *sigma2 * rho *
	      log1p(dist[currentPair] *dist[currentPair] / ( *range * *range)) /
	      mahalDist[currentPair] * jacCommon;
	    break;
	  case 3:
	    //i.e. powered exponential
	    hess[(2 * *nObs + k) * nPairs + currentPair] = -*sigma2 * rho * R_pow(dist[currentPair] / *range, *smooth) *
	      *smooth / (*range * mahalDist[currentPair]) * jacCommon;
	    hess[(3 * *nObs + k) * nPairs + currentPair] = *sigma2 * rho *
	      R_pow(dist[currentPair] / *range, *smooth) * log(dist[currentPair] / *range) /
	      mahalDist[currentPair] * jacCommon;
	    break;
	  case 4:
	    //i.e. Bessel
	    hess[(2 * *nObs + k) * nPairs + currentPair] = -*sigma2 * sill *
	      R_pow(2 * *range / dist[currentPair], *smooth) * gammafn(*smooth + 1) *
	      dist[currentPair] / (*range * *range * mahalDist[currentPair]) *
	      bessel_j(dist[currentPair] / *range, *smooth + 1) * jacCommon;
	    /* There's no closed form for the partial derivative of the
	       BesselJ function w.r.t. smooth. We use finite differences
	       for this... */
	    hess[(3 * *nObs + k) * nPairs + currentPair] = -*sigma2 * rho *
	      (M_LN2 - log(dist[currentPair] / *range) + digamma(*smooth + 1) +
	       (bessel_j(dist[currentPair] / *range, *smooth + h) -
		bessel_j(dist[currentPair] / *range, *smooth - h)) /
	       (2 * h * bessel_j(dist[currentPair] / *range, *smooth))) /
	      mahalDist[currentPair] * jacCommon;
	    break;
	  case 5:
	    //i.e. Generalized Cauchy
	    hess[(2 * *nObs + k) * nPairs + currentPair] = -*sigma2 * rho * *smooth /
	      (*range * (1 + R_pow(*range / dist[currentPair], *smooth2)) * mahalDist[currentPair]) *
	      jacCommon;
	    hess[(3 * *nObs + k) * nPairs + currentPair] = *sigma2 * log1p(R_pow(dist[currentPair] / *range, *smooth2)) *
	      rho / (*smooth2 * mahalDist[currentPair])* jacCommon;
	    hess[(4 * *nObs + k) * nPairs + currentPair] = -*sigma2 * (log1p(R_pow(dist[currentPair] / *range, *smooth2)) /
								       *smooth2 - log(dist[currentPair] / *range) /
								       (1 + R_pow(*range / dist[currentPair], *smooth2))) *
	      *smooth / *smooth2 * rho /mahalDist[currentPair] * jacCommon;
	    grad[4 * *nObs + k] += hess[(4 * *nObs + k) * nPairs + currentPair];
	    //The caugen has one more smooth parameter so
	    nCorPar = 5;
	    break;
	  }
	  grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];
	  grad[3 * *nObs + k] += hess[(3 * *nObs + k) * nPairs + currentPair];
	}
      }
    }
  }

  if (*fitmarge)
    marginalPartSmith(&nCorPar, nObs, nSite, data, frech, mahalDist, locs, scales, shapes,
		      trendlocs, trendscales, trendshapes, nloccoeff, nscalecoeff, nshapecoeff,
		      ntemploccoeff, ntempscalecoeff, ntempshapecoeff, locdsgnmat, scaledsgnmat,
		      shapedsgnmat, tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, weights,
		      hess, grad);

  return;
}

void brownresnickstderr(double *data, double *dist, int *nSite, int *nObs, double *locdsgnmat,
			int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
			double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
			int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
			double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
			double *scalecoeff, double *shapecoeff, double *temploccoeff,
			double *tempscalecoeff, double *tempshapecoeff, double *range,
			double *smooth, int *fitmarge, int *usetempcov, double *weights, double *hess,
			double *grad){
  /* This function computes the hessian of the log-pairwise
     likelihood for the Brown-Resnick model.

     Remember that this model shares the same bivariate density as the
     Smith model except that the Mahalanobis distance is modified. */

  const int nPairs = *nSite * (*nSite - 1) / 2,
    flag = usetempcov[0] + usetempcov[1] + usetempcov[2];
  int i, currentPair = -1;
  double *mahalDist, *locs, *scales, *shapes, *trendlocs, *trendscales, *trendshapes,
    *jac, *frech;

  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  trendlocs = (double *)R_alloc(*nObs, sizeof(double));
  trendscales = (double *)R_alloc(*nObs, sizeof(double));
  trendshapes = (double *)R_alloc(*nObs, sizeof(double));

  /* We have to intialized them to 0 as we don't know if
     dsgnmat2temptrend will be called */
  for (i=*nObs;i--;)
    trendlocs[i] = trendscales[i] = trendshapes[i] = 0;

  //Stage 0: Compute the covariance at each location
  brownResnick(dist, nPairs, *range, *smooth, mahalDist);

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){

    dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff, shapecoeff,
		  *nSite, *nloccoeff, *nscalecoeff, *nshapecoeff, locs, scales, shapes);

    //Stage 1: Transformation to unit Frechet
    if (flag) {
      dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, temploccoeff,
			tempscalecoeff, tempshapecoeff, *nSite, *nObs, usetempcov, *ntemploccoeff,
			*ntempscalecoeff, *ntempshapecoeff, trendlocs, trendscales, trendshapes);

      gev2frechTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs, trendscales,
		     trendshapes, jac, frech);
    }

    else
      gev2frech(data, *nObs, *nSite, locs, scales, shapes, jac, frech);

  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];

  //Stage 2: Hessian computations;
  // a- Covariance part
  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){
      currentPair++;

      int k;
      double imahal = 1 / mahalDist[currentPair], imahalSquare = imahal * imahal;

      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    hess[k * nPairs + currentPair] = hess[(*nObs + k) * nPairs + currentPair] = NA_REAL;

	    continue;
	  }

	  double ifrech1 = 1 / frech[k + i * *nObs], ifrech2 = 1 / frech[k + j * *nObs],
	    ifrech1Square = ifrech1 * ifrech1, ifrech2Square = ifrech2 * ifrech2,
	    c1 = log(frech[k + j * *nObs] * ifrech1) * imahal + 0.5 * mahalDist[currentPair],
	    c2 = mahalDist[currentPair] - c1,
	    dnormc1 = dnorm(c1, 0., 1., 0), pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0), pnormc2 = pnorm(c2, 0., 1., 1, 0),
	    //A = - pnormc1 * ifrech1 - pnormc2 * ifrech2,
	    B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	    dnormc2 * imahal * ifrech2Square,
	    C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	    dnormc1 * imahal * ifrech1Square,
	    D = c2 * dnormc1 * ifrech2 * imahalSquare * ifrech1Square +
	    c1 * dnormc2 * ifrech1 * imahalSquare * ifrech2Square,
	    dAa = - c2 * dnormc1 * ifrech1 * imahal - c1 * dnormc2 * imahal * ifrech2,
	    dBa = (c1 * c1 - 1) * dnormc2 * imahalSquare * ifrech2Square +
	    (1 + c1 * c2 ) * dnormc1 * ifrech1 * ifrech2 * imahalSquare,
	    dCa = (c2 * c2 - 1) * dnormc1 * imahalSquare * ifrech1Square +
	    (1 + c1 * c2) * dnormc2 * ifrech1 * ifrech2 * imahalSquare,
	    dDa = (c1 - c1 * c2 * c2 - 2 * c2) * dnormc1 * imahalSquare * imahal *
	    ifrech1Square * ifrech2 + (c2 - c1 * c1 *c2 - 2 * c1) * dnormc2 *
	    imahalSquare * imahal * ifrech1 * ifrech2Square,
	    jacCommon = weights[currentPair] * (dAa + (dBa * C + B * dCa + dDa) / (B*C + D));

	  hess[k * nPairs + currentPair] = -0.5 * *smooth * mahalDist[currentPair] / *range * jacCommon;
	  hess[(*nObs + k) * nPairs + currentPair] = 0.5 * log(dist[currentPair] / *range) *
	    mahalDist[currentPair] * jacCommon;

	  grad[k] += hess[k * nPairs + currentPair];
	  grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];

	}
      }
    }
  }

  if (*fitmarge){
    int start = 2;
    marginalPartSmith(&start, nObs, nSite, data, frech, mahalDist, locs, scales, shapes,
		      trendlocs, trendscales, trendshapes, nloccoeff, nscalecoeff, nshapecoeff,
		      ntemploccoeff, ntempscalecoeff, ntempshapecoeff, locdsgnmat, scaledsgnmat,
		      shapedsgnmat, tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, weights,
		      hess, grad);
  }

  return;
}


void spatgevstderr(double *data, int *nSite, int *nObs, double *locdsgnmat,
		   int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		   double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
		   int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
		   double *tempdsgnmatshape, int *ntempshapecoeff,  double *loccoeff,
		   double *scalecoeff, double *shapecoeff, double *temploccoeff,
		   double *tempscalecoeff, double *tempshapecoeff, int *usetempcov,
		   double *hess, double *grad){

  int i, flag = usetempcov[0] + usetempcov[1] + usetempcov[2];
  double *locs, *scales, *shapes, *trendlocs, *trendscales, *trendshapes;

  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  trendlocs = (double *)R_alloc(*nObs, sizeof(double));
  trendscales = (double *)R_alloc(*nObs, sizeof(double));
  trendshapes = (double *)R_alloc(*nObs, sizeof(double));

  /* We have to intialized them to 0 as we don't know if
     dsgnmat2temptrend will be called */
  for (i=*nObs;i--;)
    trendlocs[i] = trendscales[i] = trendshapes[i] = 0;

  //Stage 0. Compute the GEV parameters using the design matrix
  dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff, shapecoeff,
		*nSite, *nloccoeff, *nscalecoeff, *nshapecoeff, locs, scales, shapes);

  if (flag)
    dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, temploccoeff,
		      tempscalecoeff, tempshapecoeff, *nSite, *nObs, usetempcov, *ntemploccoeff,
		      *ntempscalecoeff, *ntempshapecoeff, trendlocs, trendscales, trendshapes);

  //Stage 1. Compute the gradient
  for (i=0;i<*nObs;i++){
    int j;
    for (j=0;j<*nSite;j++){
      int idx = *nloccoeff, k;

      if (ISNA(data[i + j * *nObs])){
	for (k=0;k<(*nloccoeff + *nscalecoeff + *nshapecoeff + *ntemploccoeff + *ntempscalecoeff + *ntempshapecoeff);k++)
	  hess[(k * *nObs + i) * *nSite + j] = NA_REAL;

	continue;
      }

      double loc = locs[j] + trendlocs[i], scale = scales[j] + trendscales[i],
	shape = shapes[j] + trendshapes[i], dataTrans = 1 + shape * (data[j * *nObs + i] - loc) / scale;

      for (k=0;k<*nloccoeff;k++){
	hess[(k * *nObs + i) * *nSite + j] = ((1 + shape) / dataTrans - R_pow(dataTrans, - 1 / shape - 1)) *
	  locdsgnmat[k * *nSite + j] / scale;
	grad[k * *nObs + i] += hess[(k * *nObs + i) * *nSite + j];
      }

      for (k=0;k<*nscalecoeff;k++){
	hess[((idx + k) * *nObs + i) * *nSite + j] =  (-1 + (1 + shape) * (data[j * *nObs + i] - loc) /
						       (scale * dataTrans) - R_pow(dataTrans, - 1 / shape - 1) * (data[j * *nObs + i] - loc) / scale) *
	  scaledsgnmat[k * *nSite + j] / scale;
	grad[(idx + k) * *nObs + i] += hess[((idx + k) * *nObs + i) * *nSite + j];
      }

      idx += *nscalecoeff;

      for (k=0;k<*nshapecoeff;k++){
	hess[((idx + k) * *nObs + i) * *nSite + j] = (log1p(dataTrans - 1) / shape - (1 + shape) * (data[j * *nObs + i] - loc) /
						      (scale * dataTrans) - R_pow(dataTrans, - 1 / shape) *
						      (log1p(dataTrans - 1) / shape - (data[j * *nObs + i] - loc) /
						       (scale * dataTrans))) * shapedsgnmat[k * *nSite + j] / shape;
	grad[(idx + k) * *nObs + i] += hess[((idx + k) * *nObs + i) * *nSite + j];
      }

      idx += *nshapecoeff;

      for (k=0;k<*ntemploccoeff;k++){
	hess[((idx + k) * *nObs + i) * *nSite + j] = ((1 + shape) / dataTrans - R_pow(dataTrans, - 1 / shape - 1)) *
	  tempdsgnmatloc[k * *nObs + i] / scale;
	grad[(idx + k) * *nObs + i] +=  hess[((idx + k) * *nObs + i) * *nSite + j];
      }

      idx += *ntemploccoeff;

      for (k=0;k<*ntempscalecoeff;k++){
	hess[((idx + k) * *nObs + i) * *nSite + j] = (-1 + (1 + shape) * (data[j * *nObs + i] - loc) /
						      (scale * dataTrans) - R_pow(dataTrans, - 1 / shape - 1) * (data[j * *nObs + i] - loc) / scale) *
	  tempdsgnmatscale[k * *nObs + i] / scale;
	grad[(idx + k) * *nObs + i] += hess[((idx + k) * *nObs + i) * *nSite + j];
      }

      idx += *ntempscalecoeff;

      for (k=0;k<*ntempshapecoeff;k++){
	hess[((idx + k) * *nObs + i) * *nSite + j] = (log1p(dataTrans - 1) / shape - (1 + shape) * (data[j * *nObs + i] - loc) /
						      (scale * dataTrans) - R_pow(dataTrans, - 1 / shape) *
						      (log1p(dataTrans - 1) / shape - (data[j * *nObs + i] - loc) /
						       (scale * dataTrans))) * tempdsgnmatshape[k * *nObs + i] / shape;
	grad[(idx + k) * *nObs + i] += hess[((idx + k) * *nObs + i) * *nSite + j];
      }
    }
  }

  return;
}

void extremaltstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		     double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		     double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
		     int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
		     double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
		     double *scalecoeff, double *shapecoeff, double *temploccoeff,
		     double *tempscalecoeff, double *tempshapecoeff, double *nugget, double *range,
		     double *smooth, double *smooth2, double *df, int *fitmarge, int *usetempcov,
		     double *weights, double *hess, double *grad){
  /* This function computes the hessian of the log-pairwise
     likelihood for the extremal t model. */

  const int nPairs = *nSite * (*nSite - 1) / 2,
    flag = usetempcov[0] + usetempcov[1] + usetempcov[2];
  int i, currentPair = -1, nCorPar = 3;
  const double idf = 1 / *df, dfPlus1 = *df + 1, eps = 0.001,
    didfidfdfPlus1_df = - (*df + 2) * idf * idf * idf,
    sill = 1 - *nugget;

  double h = 0, *jac = (double *)R_alloc(*nObs * *nSite, sizeof(double)),
    *rho = (double *)R_alloc(nPairs, sizeof(double)),
    *locs = (double *)R_alloc(*nSite, sizeof(double)),
    *scales = (double *)R_alloc(*nSite, sizeof(double)),
    *shapes = (double *)R_alloc(*nSite, sizeof(double)),
    *frech = (double *)R_alloc(*nObs * *nSite, sizeof(double)),
    *trendlocs = (double *)R_alloc(*nObs, sizeof(double)),
    *trendscales = (double *)R_alloc(*nObs, sizeof(double)),
    *trendshapes = (double *)R_alloc(*nObs, sizeof(double));


  if (*covmod == 1){
    h = sqrt(EPS) * *smooth;
    double temp = *smooth + h;
    h = temp - *smooth;
  }

  else if (*covmod == 4){
    h = sqrt(EPS) * (*smooth + 1);
    double temp = (*smooth + 1) + h;
    h = temp - (*smooth + 1);
  }

  /* We have to intialized them to 0 as we don't know if
     dsgnmat2temptrend will be called */
  for (i=*nObs;i--;)
    trendlocs[i] = trendscales[i] = trendshapes[i] = 0;

  //Stage 0: Compute the covariance at each location
  switch (*covmod){
  case 1:
    whittleMatern(dist, nPairs, *nugget, sill, *range, *smooth, rho);
    break;
  case 2:
    cauchy(dist, nPairs, *nugget, sill, *range, *smooth, rho);
    break;
  case 3:
    powerExp(dist, nPairs, *nugget, sill, *range, *smooth, rho);
    break;
  case 4:
    //Here we use 0 for dim as we don't care for the computation of
    //the hessian
    bessel(dist, nPairs, 0, *nugget, sill, *range, *smooth, rho);
    break;
  case 5:
    caugen(dist, nPairs, *nugget, sill, *range, *smooth, *smooth2, rho);
  }

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){

    dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff, shapecoeff,
		  *nSite, *nloccoeff, *nscalecoeff, *nshapecoeff, locs, scales, shapes);

    //Stage 1: Transformation to unit Frechet
    if (flag) {
      dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, temploccoeff,
			tempscalecoeff, tempshapecoeff, *nSite, *nObs, usetempcov, *ntemploccoeff,
			*ntempscalecoeff, *ntempshapecoeff, trendlocs, trendscales, trendshapes);

      gev2frechTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs, trendscales,
		     trendshapes, jac, frech);
    }

    else
      gev2frech(data, *nObs, *nSite, locs, scales, shapes, jac, frech);

  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];

  //Stage 2: Hessian computations;
  // a- Covariance part
  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){
      currentPair++;

      int k;
      double a = sqrt(dfPlus1 / (1 - rho[currentPair] * rho[currentPair])),
	//partial derivatives of a w.r.t. rho and df
	da_rho = a * rho[currentPair] / (1 - rho[currentPair] * rho[currentPair]),
	da_df = 0.5 * a / dfPlus1;

      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    hess[k * nPairs + currentPair] = hess[(*nObs + k) * nPairs + currentPair] =
	      hess[(2 * *nObs + k) * nPairs + currentPair] = hess[(3 * *nObs + k) * nPairs + currentPair] =
	      NA_REAL;

	    if (*covmod == 5)
	      hess[(4 * *nObs + k) * nPairs + currentPair] = NA_REAL;

	    continue;
	  }

	  /*-------------------------------------------------------------------------
	    We start by computing the gradient for the correlation
	    function parameters
	    --------------------------------------------------------------------------*/

	  double //some useful quantities
	    ifrech1 = 1 / frech[k + i * *nObs],
	    ifrech2 = 1 / frech[k + j * *nObs],
	    frech2_1 = R_pow(frech[k + j * *nObs] * ifrech1, idf),
	    frech1_2 = 1 / frech2_1,
	    c1 = (frech2_1 - rho[currentPair]) * a,
	    c2 = (frech1_2 - rho[currentPair]) * a,
	    //The t density and distribution evaluated at c1 and c2
	    tc1 = dt(c1, dfPlus1, 0),
	    tc2 = dt(c2, dfPlus1, 0),
	    Tc1 = pt(c1, dfPlus1, 1, 0),
	    Tc2 = pt(c2, dfPlus1, 1, 0),
	    //the next following four variables are the first and second
	    //derivative of the t density evaluated at c1 and c2
	    dertc1 = -(dfPlus1 + 1) * c1 / (dfPlus1 + c1 * c1) * tc1,
	    dertc2 = -(dfPlus1 + 1) * c2 / (dfPlus1 + c2 * c2) * tc2,
	    der2tc1 = dertc1 * dertc1 / tc1 + dertc1 / c1 + 2 * dertc1 * dertc1 /
	    (tc1 * (dfPlus1 + 1)),
	    der2tc2 = dertc2 * dertc2 / tc2 + dertc2 / c2 + 2 * dertc2 * dertc2 /
	    (tc2 * (dfPlus1 + 1)),
	    //partial derivatives of c1, c2 w.r.t. rho
	    dc1_rho = -a + c1 / a * da_rho,
	    dc2_rho = -a + c2 / a * da_rho,
	    //A = -Tc1 * ifrech1 - Tc2 * ifrech2,
	    dA_rho = -ifrech1 * dc1_rho * tc1 - ifrech2 * dc2_rho * tc2,
	    B = ifrech1 * ifrech1 * Tc1 + ifrech1 * ifrech1 * idf * frech2_1 * tc1 * a -
	    ifrech1 * ifrech2 * idf * frech1_2 * tc2 * a,
	    dB_rho = ifrech1 * ifrech1 * dc1_rho * tc1 +
	    ifrech1 * ifrech1 * idf * frech2_1 * dc1_rho * dertc1 * a +
	    ifrech1 * ifrech1 * idf * frech2_1 * tc1 * da_rho -
	    ifrech1 * ifrech2 * idf * frech1_2 * dc2_rho * dertc2 * a -
	    ifrech1 * ifrech2 * idf * frech1_2 * tc2 * da_rho,
	    C = ifrech2 * ifrech2 * Tc2 + ifrech2 * ifrech2 * idf * frech1_2 * tc2 * a -
	    ifrech1 * ifrech2 * idf * frech2_1 * tc1 * a,
	    dC_rho = ifrech2 * ifrech2 * dc2_rho * tc2 +
	    ifrech2 * ifrech2 * idf * frech1_2 * dc2_rho * dertc2 * a +
	    ifrech2 * ifrech2 * idf * frech1_2 * tc2 * da_rho -
	    ifrech1 * ifrech2 * idf * frech2_1 * dc1_rho * dertc1 * a -
	    ifrech1 * ifrech2 * idf * frech2_1 * tc1 * da_rho,
	    D = ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * dfPlus1 * a * tc1 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * dfPlus1 * a * tc2 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a * dertc1 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a * dertc2,
	    dD_rho = ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * dfPlus1 * da_rho * tc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * dfPlus1 * a * dc1_rho * dertc1 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * dfPlus1 * da_rho * tc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * dfPlus1 * a * dc2_rho * dertc2 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * 2 * da_rho * a * dertc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a * dc1_rho * der2tc1 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * 2 * da_rho * a * dertc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a * dc2_rho * der2tc2,
	    iBCplusD = 1 / (B * C + D),
	    jacCommonRho = weights[currentPair] * (dA_rho + (dB_rho * C + B * dC_rho + dD_rho) * iBCplusD);

	  hess[k * nPairs + currentPair] = -rho[currentPair] / sill * jacCommonRho;
	  grad[k] += hess[k * nPairs + currentPair];

	  switch (*covmod){
	  case 1:
	    //i.e. Whittle-Matern
	    hess[(*nObs + k) * nPairs + currentPair] = rho[currentPair] *
	      (-2 * *smooth / *range + dist[currentPair] *
	       bessel_k(dist[currentPair] / *range, *smooth + 1, 1) /
	       (bessel_k(dist[currentPair] / *range, *smooth, 1) * *range * *range)) *
	      jacCommonRho;
	    /* There's no closed form for the partial derivative of the
	       BesselK function w.r.t. smooth. We use finite differences
	       for this... */
	    hess[(2 * *nObs + k) * nPairs + currentPair] = rho[currentPair] *
	      (-M_LN2 - digamma(*smooth) + log(dist[currentPair] / *range) +
	       (bessel_k(dist[currentPair] / *range, *smooth + h, 1) -
		bessel_k(dist[currentPair] / *range, *smooth - h, 1)) /
	       (2 * h * bessel_k(dist[currentPair] / *range, *smooth, 1))) *
	      jacCommonRho;
	    break;
	  case 2:
	    //i.e. cauchy
	    hess[(*nObs + k) * nPairs + currentPair] = 2 * dist[currentPair] * dist[currentPair] *
	      sill * *smooth / (*range * *range * *range) *
	      R_pow(1 + dist[currentPair] * dist[currentPair] / (*range * *range),
		    - *smooth - 1) * jacCommonRho;
	    hess[(2 * *nObs + k) * nPairs + currentPair] = -rho[currentPair] *
	      log1p(dist[currentPair] *dist[currentPair] / ( *range * *range)) *
	      jacCommonRho;
	    break;
	  case 3:
	    //i.e. powered exponential
	    hess[(*nObs + k) * nPairs + currentPair] = rho[currentPair] * *smooth / *range *
	      R_pow(dist[currentPair] / *range, *smooth) * jacCommonRho;
	    hess[(2 * *nObs + k) * nPairs + currentPair] = -rho[currentPair] *
	      R_pow(dist[currentPair] / *range, *smooth) * log(dist[currentPair] / *range) *
	      jacCommonRho;
	    break;
	  case 4:
	    //i.e. Bessel
	    hess[(*nObs + k) * nPairs + currentPair] = sill *
	      R_pow(2 * *range / dist[currentPair], *smooth) *
	      dist[currentPair] / (*range * *range) * gammafn(*smooth + 1) *
	      bessel_j(dist[currentPair] / *range, *smooth + 1) * jacCommonRho;
	    /* There's no closed form for the partial derivative of the
	       BesselJ function w.r.t. smooth. We use finite differences
	       for this... */
	    hess[(2 * *nObs + k) * nPairs + currentPair] = rho[currentPair] *
	      (M_LN2 - log(dist[currentPair] / *range) + digamma(*smooth + 1) +
	       (bessel_j(dist[currentPair] / *range, *smooth + h) -
		bessel_j(dist[currentPair] / *range, *smooth - h)) /
	       (2 * h * bessel_j(dist[currentPair] / *range, *smooth))) *
	      jacCommonRho;
	    break;
	  case 5:
	    //i.e. Generalized Cauchy
	    hess[(*nObs + k) * nPairs + currentPair] = rho[currentPair] * *smooth / *range /
	      (1 + R_pow(*range / dist[currentPair], *smooth2)) *
	      jacCommonRho;
	    hess[(2 * *nObs + k) * nPairs + currentPair] = -log1p(R_pow(dist[currentPair] / *range, *smooth2)) *
	      rho[currentPair] / *smooth2 *jacCommonRho;
	    hess[(3 * *nObs + k) * nPairs + currentPair] = (log1p(R_pow(dist[currentPair] / *range, *smooth2)) /
							    *smooth2 - log(dist[currentPair] / *range) /
							    (1 + R_pow(*range / dist[currentPair], *smooth2))) *
	      *smooth / *smooth2 * rho[currentPair] * jacCommonRho;
	    grad[3 * *nObs + k] += hess[(3 * *nObs + k) * nPairs + currentPair];
	    //The caugen has 2 smooth paramaters so add + 1 to nCorPar.
	    nCorPar = 4;
	    break;
	  }

	  grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];
	  grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];



	  /*-------------------------------------------------------------------------
	    And now compute the gradient for the degree of freedom
	    --------------------------------------------------------------------------*/
	  double //some useful quantities
	    //partial derivatives of frech2_1, frech1_2 w.r.t. df
	    dfrech21_df = - frech2_1 * log(frech2_1) * idf,
	    dfrech12_df = - frech1_2 * log(frech1_2) * idf,
	    //partial derivatives of c1, c2 w.r.t. df
	    dc1_df = a * dfrech21_df + (frech2_1 - rho[currentPair]) * da_df,
	    dc2_df = a * dfrech12_df + (frech1_2 - rho[currentPair]) * da_df,
	    /* Here we have to be aware as c1 is a function of df and
	       the t distribution/density depend also on df !!! Hence the
	       generic derivation is of the form :

	       \frac{\partial f(c1,df)}{\partial df} = D_1(f)(c1,df)
	       \frac{\partial c1}{\partial df} + D_2(f)(c1,df),

	       where D_i(f) denotes the derivative w.r.t. to the $i$-th
	       variable of the function $f$. */

	    /* No closed form exists for D_2[T](c1,df) so we use finite
	       differences for this... */
	    dTc1_df = dc1_df * tc1 + (pt(c1, dfPlus1 + eps, 1, 0) - Tc1) / eps,
	    dTc2_df = dc2_df * tc2 + (pt(c2, dfPlus1 + eps, 1, 0) - Tc2) / eps,
	    /* Below is the derivative of the t density w.r.t. df and
	       evaluated at c1 and c2 i.e. c1/c2 aren't differentiated
	       w.r.t. df!!! */
	    D_2_tc1 = (digamma(0.5 * (dfPlus1 + 1)) - 1 / dfPlus1 -
		       digamma(0.5 * dfPlus1) - log1p(c1 * c1 / dfPlus1) +
		       (dfPlus1 + 1) * c1 * c1 / (dfPlus1 * (dfPlus1 + c1 * c1))) *
	    0.5 * tc1,
	    D_2_tc2 = (digamma(0.5 * (dfPlus1 + 1)) - 1 / dfPlus1 -
		       digamma(0.5 * dfPlus1) - log1p(c2 * c2 / dfPlus1) +
		       (dfPlus1 + 1) * c2 * c2 / (dfPlus1 * (dfPlus1 + c2 * c2))) *
	    0.5 * tc2,
	    /* Below are the total derivative of tc1 w.r.t. df i.e. both
	       c1 and the t density are derived w.r.t. df */
	    dtc1_df = dc1_df * dertc1 + D_2_tc1,
	    dtc2_df = dc2_df * dertc2 + D_2_tc2,
	    /* Below are the total derivative of dertc1 and dertc2
	       w.r.t. df */
	    ddertc1_df = dertc1 / (dfPlus1 + 1) + dertc1 / tc1 * dtc1_df + dertc1 / c1 * dc1_df -
	    dertc1 / (dfPlus1 + c1 * c1) * (1 + 2 * dc1_df * c1),
	    ddertc2_df = dertc2 / (dfPlus1 + 1) + dertc2 / tc2 * dtc2_df + dertc2 / c2 * dc2_df -
	    dertc2 / (dfPlus1 + c2 * c2) * (1 + 2 * dc2_df * c2),
	    //A = -Tc1 * ifrech1 - Tc2 * ifrech2,
	    dA_df = -ifrech1 * dTc1_df - ifrech2 * dTc2_df,
	    //B = ifrech1 * ifrech1 * Tc1 + ifrech1 * ifrech1 * idf * frech2_1 * tc1 * a -
	    //ifrech1 * ifrech2 * idf * frech1_2 * tc2 * a,
	    dB_df = ifrech1 * ifrech1 * dTc1_df -
	    ifrech1 * ifrech1 * idf * idf * frech2_1 * tc1 * a +
	    ifrech1 * ifrech1 * idf * dfrech21_df * tc1 * a +
	    ifrech1 * ifrech1 * idf * frech2_1 * dtc1_df * a +
	    ifrech1 * ifrech1 * idf * frech2_1 * tc1 * da_df +
	    ifrech1 * ifrech2 * idf * idf * frech1_2 * tc2 * a -
	    ifrech1 * ifrech2 * idf * dfrech12_df * tc2 * a -
	    ifrech1 * ifrech2 * idf * frech1_2 * dtc2_df * a -
	    ifrech1 * ifrech2 * idf * frech1_2 * tc2 * da_df,
	    //C = ifrech2 * ifrech2 * Tc2 + ifrech2 * ifrech2 * idf * frech1_2 * tc2 * a -
	    //ifrech1 * ifrech2 * idf * frech2_1 * tc1 * a,
	    dC_df = ifrech2 * ifrech2 * dTc2_df -
	    ifrech2 * ifrech2 * idf * idf * frech1_2 * tc2 * a +
	    ifrech2 * ifrech2 * idf * dfrech12_df * tc2 * a +
	    ifrech2 * ifrech2 * idf * frech1_2 * dtc2_df * a +
	    ifrech2 * ifrech2 * idf * frech1_2 * tc2 * da_df +
	    ifrech1 * ifrech2 * idf * idf * frech2_1 * tc1 * a -
	    ifrech1 * ifrech2 * idf * dfrech21_df * tc1 * a -
	    ifrech1 * ifrech2 * idf * frech2_1 * dtc1_df * a -
	    ifrech1 * ifrech2 * idf * frech2_1 * tc1 * da_df,
	    //D = ifrech1 * ifrech1 * ifrech2 * idf * idf * dfPlus1 * frech2_1 * a * tc1 +
	    //ifrech1 * ifrech2 * ifrech2 * idf * idf * dfPlus1 * frech1_2 * a * tc2 +
	    //ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a * dertc1 +
	    //ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a * dertc2,
	    dD_df = ifrech1 * ifrech1 * ifrech2 * didfidfdfPlus1_df * frech2_1 * a * tc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * dfPlus1 * dfrech21_df * a * tc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * dfPlus1 * frech2_1 * da_df * tc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * dfPlus1 * frech2_1 * a * dtc1_df +
	    ifrech1 * ifrech2 * ifrech2 * didfidfdfPlus1_df * frech1_2 * a * tc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * dfPlus1 * dfrech12_df * a * tc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * dfPlus1 * frech1_2 * da_df * tc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * dfPlus1 * frech1_2 * a * dtc2_df -
	    2 * ifrech1 * ifrech1 * ifrech2 * idf * idf * idf * frech2_1 * frech2_1 * a * a *
	    dertc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * 2 * dfrech21_df * frech2_1 * a * a * dertc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * 2 * da_df * a * dertc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a * ddertc1_df -
	    2 * ifrech1 * ifrech2 * ifrech2 * idf * idf * idf * frech1_2 * frech1_2 * a * a *
	    dertc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * 2 * dfrech12_df * frech1_2 * a * a * dertc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * 2 * da_df * a * dertc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a * ddertc2_df;

	  hess[(nCorPar * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	    (dA_df + (dB_df * C + B * dC_df + dD_df) * iBCplusD);
	  grad[nCorPar * *nObs + k] += hess[(nCorPar * *nObs + k) * nPairs + currentPair];

	}
      }
    }
  }

  if (*fitmarge){
    nCorPar++;
    marginalPartExtremalt(&nCorPar, nObs, nSite, data, frech, df, rho, locs, scales, shapes,
			  trendlocs, trendscales, trendshapes, nloccoeff, nscalecoeff,
			  nshapecoeff, ntemploccoeff, ntempscalecoeff, ntempshapecoeff,
			  locdsgnmat, scaledsgnmat, shapedsgnmat, tempdsgnmatloc,
			  tempdsgnmatscale, tempdsgnmatshape, weights, hess, grad);
  }

  return;
}
