#include "header.h"

void copula(int *copula, int *covmod, double *dist, double *data, int *nSite, int *nObs,
	    int *dim, int *fitmarge, double *locdsgnmat,  double *locpenmat, int *nloccoeff,
	    int *npparloc, double *locpenalty, double *scaledsgnmat, double *scalepenmat,
	    int *nscalecoeff, int *npparscale, double *scalepenalty,
	    double *shapedsgnmat, double *shapepenmat, int *nshapecoeff, int *npparshape,
	    double *shapepenalty, int *usetempcov, double *tempdsgnmatloc,
	    double *temppenmatloc, int *ntempcoeffloc, int *nppartempcoeffloc,
	    double *temppenaltyloc, double *tempdsgnmatscale, double *temppenmatscale,
	    int *ntempcoeffscale, int *nppartempcoeffscale, double *temppenaltyscale,
	    double *tempdsgnmatshape, double *temppenmatshape, int *ntempcoeffshape,
	    int *nppartempcoeffshape, double *temppenaltyshape, double *loccoeff,
	    double *scalecoeff, double *shapecoeff, double *tempcoeffloc,
	    double *tempcoeffscale, double *tempcoeffshape, double *DoF, double *nugget, double *range,
	    double *smooth, double *smooth2, double *dns){

  int flag = usetempcov[0] + usetempcov[1] + usetempcov[2],
    nPairs = *nSite * (*nSite + 1) / 2;

  if (*nugget >= 1){
    *dns = (*nugget * *nugget) * MINF;
    return;
  }

  double *trendlocs = malloc(*nObs * sizeof(double)),
    *trendscales = malloc(*nObs * sizeof(double)),
    *trendshapes = malloc(*nObs * sizeof(double)),
    *logdens = malloc(*nSite * *nObs * sizeof(double)),
    *covariances = malloc(nPairs * sizeof(double)),
    *covMat = malloc(*nSite * *nSite * sizeof(double)),
    *locs = malloc(*nSite * sizeof(double)),
    *scales = malloc(*nSite * sizeof(double)),
    *shapes = malloc(*nSite * sizeof(double)),
    *unif = malloc(*nSite * *nObs * sizeof(double)),
    sill = 1 - *nugget;

  for (int i=0;i<(*nSite * *nSite);i++)
    covMat[i] = 0;

  //Stage 1: Compute the covariance at each location
  switch (*covmod){
  case 1:
    *dns = whittleMatern(dist, nPairs, *nugget, sill, *range, *smooth, covariances);
    break;
  case 2:
    *dns = cauchy(dist, nPairs, *nugget, sill, *range, *smooth, covariances);
    break;
  case 3:
    *dns = powerExp(dist, nPairs, *nugget, sill, *range, *smooth, covariances);
    break;
  case 4:
    *dns = bessel(dist, nPairs, *dim, *nugget, sill, *range, *smooth, covariances);
    break;
  case 5:
    *dns = caugen(dist, nPairs, *nugget, sill, *range, *smooth, *smooth2, covariances);
    break;
  }

  if (*dns != 0.0)
    return;

  /* We need to fill in the upper triangular part of covMatChol with
     covariances */
  {
    int current=-1;
    for (int i=0;i<*nSite;i++)
      for (int j=i;j<*nSite;j++){
	current++;
	covMat[i + j * *nSite] = covariances[current];
      }
  }

  if (*fitmarge){
    //Stage 2: Compute the GEV parameters using the design matrix
    *dns = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff,
			 scalecoeff, shapecoeff, *nSite, *nloccoeff, *nscalecoeff,
			 *nshapecoeff, locs, scales, shapes);

    if (flag){
      dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape,
			tempcoeffloc, tempcoeffscale, tempcoeffshape, *nSite,
			*nObs, usetempcov, *ntempcoeffloc, *ntempcoeffscale,
			*ntempcoeffshape, trendlocs, trendscales, trendshapes);

      for (int i=0;i<*nSite;i++)
	for (int j=0;j<*nObs;j++)
	  if (((scales[i] + trendscales[j]) <= 0) ||
	      ((shapes[i] + trendshapes[j]) <= -1)){
	    *dns = MINF;
	    return;
	  }
    }

    else if (*dns != 0.0)
      return;
  }

  else
    for (int i=0;i<*nSite;i++)
      locs[i] = scales[i] = shapes[i] = 1;

  //Stage 3: Transformation to unit Frechet
  if (flag)
    *dns = gev2unifTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs,
			 trendscales, trendshapes, unif, logdens);

  else
    *dns = gev2unif(data, *nObs, *nSite, locs, scales, shapes, unif, logdens);

  if (*dns != 0.0)
    return;

  //Stage 4: Computation of the log-likelihood

  // First the "dependence" log-likelihood
  if (*copula == 1){
    for (int i=0;i<(*nSite * *nObs);i++)
      unif[i] = qnorm(unif[i], 0, 1, 1, 0);

    *dns = gaussianCopula(unif, 1, covMat, *nObs, *nSite);
  }

  else if (*copula == 2){

    if (*DoF <= 0){
      *dns = (1 - *DoF) * (1 - *DoF) * MINF;
      return;
    }

    if (*DoF > 35){
      *dns = (*DoF - 34) * (*DoF - 34) * MINF;
      return;
    }

    for (int i=0;i<(*nSite * *nObs);i++)
      unif[i] = qt(unif[i], *DoF, 1, 0);

    *dns = studentCopula(unif, *DoF, covMat, *nObs, *nSite);
  }

  if (*dns == MINF)
    return;

  /* Second the "marginal" log-likelihood i.e., log of the GEV
  densities */
  for (int i=0;i<(*nSite * *nObs);i++)
    *dns += logdens[i];

  //Stage 5: Removing the penalizing terms (if any)
  // 1- For the location parameter
  if (*locpenalty > 0)
    *dns -= penalization(locpenmat, loccoeff, *locpenalty, *nloccoeff,
			 *npparloc);

  // 2- For the scale parameter
  if (*scalepenalty > 0)
    *dns -= penalization(scalepenmat, scalecoeff, *scalepenalty, *nscalecoeff,
			 *npparscale);

  // 3- For the shape parameter
  if (*shapepenalty > 0)
    *dns -= penalization(shapepenmat, shapecoeff, *shapepenalty, *nshapecoeff,
			 *npparshape);

  // 4- Doing the same thing for the temporal component
  if (*temppenaltyloc > 0)
    *dns -= penalization(temppenmatloc, tempcoeffloc, *temppenaltyloc,
			 *ntempcoeffloc, *nppartempcoeffloc);

  if (*temppenaltyscale > 0)
    *dns -= penalization(temppenmatscale, tempcoeffscale, *temppenaltyscale,
			 *ntempcoeffscale, *nppartempcoeffscale);

  if (*temppenaltyshape > 0)
    *dns -= penalization(temppenmatshape, tempcoeffshape, *temppenaltyshape,
			 *ntempcoeffshape, *nppartempcoeffshape);

  if (!R_FINITE(*dns))
    *dns = MINF;


  free(trendlocs); free(trendscales); free(trendshapes); free(logdens); free(covariances);
  free(covMat); free(locs); free(scales); free(shapes); free(unif);
  return;
}

double gaussianCopula(double *data, double sd, double *covMat, int nObs,
		      int nSite){
  // This function computes the log-likelihood for the Gaussian copula

  int info = 0, oneInt = 1;
  double ans = 0, logDet = 0, one = 1;

  // Cholesky decomposition
  F77_CALL(dpotrf)("U", &nSite, covMat, &nSite, &info FCONE);

  if (info != 0)
    return MINF;

  // Compute the log of the determinant from this Cholesky decomp
  for (int i=0;i<nSite;i++)
    logDet += log(covMat[i * (nSite + 1)]);

  logDet *= 2;

  ans = - 0.5 * nObs * (logDet + nSite * log(M_2PI));

  // Compupte the log of the multivariate Gaussian density
  double *dummy = malloc(nSite * sizeof(double));
  for (int i=0;i<nObs;i++){
    for (int j=0;j<nSite;j++)
      dummy[j] = data[i + j * nObs];

    F77_CALL(dtrsm)("L", "U", "T", "N", &nSite, &oneInt, &one, covMat,
		    &nSite, dummy, &nSite
		    FCONE FCONE FCONE FCONE);

    for (int j=0;j<nSite;j++)
      ans -= 0.5 * dummy[j] * dummy[j];
  }

  //Jacobian part of the gaussian copula
  for (int i=0;i<(nSite * nObs);i++)
    ans -= dnorm(data[i], 0, sd, 1);

  free(dummy);
  return ans;
}


double studentCopula(double *data, double DoF, double *covMat, int nObs,
		     int nSite){
  // This function computes the log-likelihood for the Student copula

  int info = 0, oneInt = 1;
  double logDet = 0, one = 1, iDoF = 1 / DoF;

  // Cholesky decomposition
  F77_CALL(dpotrf)("U", &nSite, covMat, &nSite, &info FCONE);

  if (info != 0)
    return MINF;

  // Compute the log of the determinant from this Cholesky decomp
  for (int i=0;i<nSite;i++)
    logDet += log(covMat[i * (nSite + 1)]);

  logDet *= 2;

  // Compupte the log of the multivariate Student density
  double *dummy = malloc(nSite * sizeof(double)),
    ans = 0;

  for (int i=0;i<nObs;i++){
    double dummy2 = 0;

    for (int j=0;j<nSite;j++)
      dummy[j] = data[i + j * nObs];

    F77_CALL(dtrsm)("L", "U", "T", "N", &nSite, &oneInt, &one, covMat,
		    &nSite, dummy, &nSite
		    FCONE FCONE FCONE FCONE);

    for (int j=0;j<nSite;j++)
      dummy2 += dummy[j] * dummy[j];

    ans += log1p(dummy2 * iDoF);
  }

  ans = nObs * (lgammafn(0.5 * (DoF + nSite)) - lgammafn(0.5 * DoF) -
		  0.5 * nSite * log(M_PI * DoF) - 0.5 * logDet) - 0.5 *
    (DoF + nSite) * ans;

  //Jacobian part of the Student copula
  for (int i=0;i<(nSite * nObs);i++)
    ans -= dt(data[i], DoF, 1);

  free(dummy);
  return ans;
}
