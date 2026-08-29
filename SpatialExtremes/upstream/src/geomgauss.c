#include "header.h"

void geomgaussfull(int *covmod, double *data, double *dist, int *nSite,
		   int *nObs, int *dim, int *weighted, double *weights,
		   double *locs, double *scales, double *shapes,
		   double *sigma2, double *sigma2Bound, double *nugget,
		   double *range, double *smooth, double *smooth2,
		   int *fitmarge, double *dns){
  //This is the geometric gaussian model. It computes the pairwise
  //log-likelihood

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

   double *jac = malloc(*nSite * *nObs * sizeof(double)),
    *rho = malloc(nPairs * sizeof(double)),
     *frech = malloc(*nSite * *nObs * sizeof(double));

  //Stage 0: Compute the covariance at each location
  *dns = geomCovariance(dist, nPairs, *dim, *covmod, *sigma2, *sigma2Bound,
			*nugget, *range, *smooth, *smooth2, rho);

  if (*dns != 0.0){
    free(jac); free(rho); free(frech);
    return;
  }

  //Stage 1: Transformation to unit Frechet
  if (*fitmarge){
    *dns = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);

    if (*dns != 0.0){
      free(jac); free(rho); free(frech);
      return;
    }

    if (*weighted)
      *dns = wlpliksmith(frech, rho, jac, *nObs, *nSite, weights);

    else
      *dns = lpliksmith(frech, rho, jac, *nObs, *nSite);
  }

  else {
    for (int i=(*nSite * *nObs);i--;)
      jac[i] = 0;

    if (*weighted)
      *dns = wlpliksmith(data, rho, jac, *nObs, *nSite, weights);

    else
      *dns = lpliksmith(data, rho, jac, *nObs, *nSite);
  }

  if (!R_FINITE(*dns))
    *dns = MINF;

  free(jac); free(rho); free(frech);
  return;

}

void geomgaussdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
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
		      double *sigma2, double *sigma2Bound, double *nugget, double *range,
		      double *smooth, double *smooth2, double *dns){
  //This is the geometric gaussian model
  //The GEV parameters are defined using a polynomial response surface

  const int nPairs = *nSite * (*nSite - 1) / 2;
  int flag = usetempcov[0] + usetempcov[1] + usetempcov[2];

  if (*nugget >= 1){
    *dns = *nugget * *nugget * MINF;
    return;
  }

  //Stage 1: Compute the covariance at each location
  double *rho = malloc(nPairs * sizeof(double));
  *dns = geomCovariance(dist, nPairs, *dim, *covmod, *sigma2, *sigma2Bound, *nugget, *range,
			*smooth, *smooth2, rho);

  if (*dns != 0.0){
    free(rho);
    return;
  }

  double *locs = malloc(*nSite * sizeof(double)),
    *scales = malloc(*nSite * sizeof(double)),
    *shapes = malloc(*nSite * sizeof(double)),
    *trendlocs = malloc(*nObs * sizeof(double)),
    *trendscales = malloc(*nObs * sizeof(double)),
    *trendshapes = malloc(*nObs * sizeof(double));

  //Stage 2: Compute the GEV parameters using the design matrix
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
	  free(locs); free(scales); free(shapes); free(trendlocs); free(trendscales);
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
    free(trendshapes); free(jac); free(frech);
    return;
  }

  if (*weighted)
    *dns = wlpliksmith(frech, rho, jac, *nObs, *nSite, weights);

  else
    *dns = lpliksmith(frech, rho, jac, *nObs, *nSite);

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

  free(rho); free(locs); free(scales); free(shapes); free(trendlocs); free(trendscales);
  free(trendshapes); free(jac); free(frech);
  return;

}
