#include "header.h"

void schlatherfull(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		   int *dim, int *weighted, double *weights, double *locs,
		   double *scales, double *shapes, double *nugget, double *range,
		   double *smooth, double *smooth2, int *fitmarge,double *dns){
  //This is the schlather's model. It's a wrapper to several
  //sub-functions. It's named xxxfull as it either assume that the
  //margins are unit Frechet, or the GEV parameters are estimated at
  //each locations.

  const int nPairs = *nSite * (*nSite - 1) / 2;

  //Some preliminary steps: Valid points?
  if (*fitmarge){
    for (int i=0;i<*nSite;i++){
      if ((scales[i] <= 0) || (shapes[i] <= -1)){
	*dns = MINF;
	return;
      }
    }
  }

  if (*nugget >= 1){
    *dns = *nugget * *nugget * MINF;
    return;
  }

  double *rho = malloc(nPairs * sizeof(double));

  //Stage 0: Compute the covariance at each location
  switch (*covmod){
  case 1:
    *dns = whittleMatern(dist, nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 2:
    *dns = cauchy(dist, nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 3:
    *dns = powerExp(dist, nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 4:
    *dns = bessel(dist, nPairs, *dim, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 5:
    *dns = caugen(dist, nPairs, *nugget, 1 - *nugget, *range, *smooth, *smooth2, rho);
    break;
  }

  if (*dns != 0.0){
    free(rho);
    return;
  }

  //Stage 1: Transformation to unit Frechet
  double *jac = malloc(*nSite * *nObs * sizeof(double)),
    *frech = malloc(*nSite * *nObs * sizeof(double));

  if (*fitmarge){
    *dns = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);

    if (*dns != 0.0){
      free(rho); free(jac); free(frech);
      return;
    }

    if (*weighted)
      *dns = wlplikschlather(frech, rho, jac, *nObs, *nSite, weights);

    else
      *dns = lplikschlather(frech, rho, jac, *nObs, *nSite);
  }

  else {
    for (int i=0;i<(*nSite * *nObs);i++)
      jac[i] = 0;

    if (*weighted)
      *dns = wlplikschlather(data, rho, jac, *nObs, *nSite, weights);

    else
      *dns = lplikschlather(data, rho, jac, *nObs, *nSite);
  }

  free(rho); free(jac); free(frech);
  return;

}

void schlatherdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
		      int *weighted, double *weights, double *locdsgnmat, double *locpenmat,
		      int *nloccoeff, int *npparloc, double *locpenalty, double *scaledsgnmat,
		      double *scalepenmat, int *nscalecoeff, int *npparscale,
		      double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		      int *nshapecoeff, int *npparshape, double *shapepenalty, int *usetempcov,
		      double *tempdsgnmatloc, double *temppenmatloc, int *ntempcoeffloc,
		      int *nppartempcoeffloc, double *temppenaltyloc, double *tempdsgnmatscale,
		      double *temppenmatscale, int *ntempcoeffscale, int *nppartempcoeffscale,
		      double *temppenaltyscale, double *tempdsgnmatshape, double *temppenmatshape,
		      int *ntempcoeffshape, int *nppartempcoeffshape, double *temppenaltyshape,
		      double *loccoeff, double *scalecoeff, double *shapecoeff,
		      double *tempcoeffloc, double *tempcoeffscale, double *tempcoeffshape,
		      double *nugget, double *range, double *smooth, double *smooth2, double *dns){
  //This is the Schlather's model.
  //The GEV parameters are defined using a polynomial response surface

  const int nPairs = *nSite * (*nSite - 1) / 2;
  int flag = usetempcov[0] + usetempcov[1] + usetempcov[2];

  if (*nugget >= 1){
    *dns = *nugget * *nugget * MINF;
    return;
  }

  double *rho = malloc(nPairs * sizeof(double));

  //Stage 1: Compute the covariance at each location
  switch (*covmod){
  case 1:
    *dns = whittleMatern(dist, nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 2:
    *dns = cauchy(dist, nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 3:
    *dns = powerExp(dist, nPairs, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 4:
    *dns = bessel(dist, nPairs, *dim, *nugget, 1 - *nugget, *range, *smooth, rho);
    break;
  case 5:
    *dns = caugen(dist, nPairs, *nugget, 1 - *nugget, *range, *smooth, *smooth2, rho);
    break;
  }

  if (*dns != 0.0){
    free(rho);
    return;
  }

  //Stage 2: Compute the GEV parameters using the design matrix
  double *locs = malloc(*nSite * sizeof(double)),
    *scales = malloc(*nSite * sizeof(double)),
    *shapes = malloc(*nSite * sizeof(double)),
    *trendlocs = malloc(*nObs * sizeof(double)),
    *trendscales = malloc(*nObs * sizeof(double)),
    *trendshapes = malloc(*nObs * sizeof(double));

  *dns = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff, shapecoeff,
		       *nSite, *nloccoeff, *nscalecoeff, *nshapecoeff, locs, scales, shapes);

  if (flag){
    dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, tempcoeffloc,
		      tempcoeffscale, tempcoeffshape, *nSite, *nObs, usetempcov, *ntempcoeffloc,
		      *ntempcoeffscale, *ntempcoeffshape, trendlocs, trendscales, trendshapes);

    for (int i=0;i<*nSite;i++)
      for (int j=0;j<*nObs;j++)
	if (((scales[i] + trendscales[j]) <= 0) || ((shapes[i] + trendshapes[j]) <= -1)){
	  *dns = MINF;
	  free(rho); free(locs); free(scales); free(shapes); free(trendlocs); free(trendscales);
	  free(trendshapes);
	  return;
	}
  }

  else if (*dns != 0.0){
    free(rho); free(locs); free(scales); free(shapes); free(trendlocs); free(trendscales);
    free(trendshapes);
    return;
  }

  //Stage 3: Transformation to unit Frechet
   double *jac = malloc(*nObs * *nSite * sizeof(double)),
     *frech = malloc(*nObs * *nSite * sizeof(double));

  if (flag)
    *dns = gev2frechTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs, trendscales,
			  trendshapes, jac, frech);

  else
    *dns = gev2frech(data, *nObs, *nSite, locs, scales, shapes, jac, frech);

  if (*dns != 0.0){
    free(rho); free(locs); free(scales); free(shapes); free(trendlocs); free(trendscales);
    free(trendshapes); free(frech); free(jac);
    return;
  }

  if (*weighted)
    *dns = wlplikschlather(frech, rho, jac, *nObs, *nSite, weights);

  else
    *dns = lplikschlather(frech, rho, jac, *nObs, *nSite);

  //Stage 5: Removing the penalizing terms (if any)
  // 1- For the location parameter
  if (*locpenalty > 0)
    *dns -= penalization(locpenmat, loccoeff, *locpenalty, *nloccoeff, *npparloc);

  // 2- For the scale parameter
  if (*scalepenalty > 0)
    *dns -= penalization(scalepenmat, scalecoeff, *scalepenalty, *nscalecoeff, *npparscale);

  // 3- For the shape parameter
  if (*shapepenalty > 0)
    *dns -= penalization(shapepenmat, shapecoeff, *shapepenalty, *nshapecoeff, *npparshape);

  // 4- Doing the same thing for the temporal component
  if (*temppenaltyloc > 0)
    *dns -= penalization(temppenmatloc, tempcoeffloc, *temppenaltyloc, *ntempcoeffloc,
			 *nppartempcoeffloc);

  if (*temppenaltyscale > 0)
    *dns -= penalization(temppenmatscale, tempcoeffscale, *temppenaltyscale, *ntempcoeffscale,
			 *nppartempcoeffscale);

  if (*temppenaltyshape > 0)
    *dns -= penalization(temppenmatshape, tempcoeffshape, *temppenaltyshape, *ntempcoeffshape,
			 *nppartempcoeffshape);

  free(locs); free(scales); free(shapes); free(trendlocs); free(trendscales);
  free(trendshapes); free(frech); free(jac); free(rho);
  return;
}
