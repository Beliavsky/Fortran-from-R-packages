#include "header.h"

void spatgevlik(double *data, double *covariables, int *nSite, int *nObs,
		double *locdsgnmat, double *locpenmat, int *nloccoeff,
		int *npparloc, double *locpenalty, double *scaledsgnmat,
		double *scalepenmat, int *nscalecoeff, int *npparscale,
		double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		int *nshapecoeff, int *npparshape, double *shapepenalty,
		int *usetempcov, double *tempdsgnmatloc, double *temppenmatloc,
		int *ntempcoeffloc, int *nppartempcoeffloc, double *temppenaltyloc,
		double *tempdsgnmatscale, double *temppenmatscale, int *ntempcoeffscale,
		int *nppartempcoeffscale, double *temppenaltyscale, double *tempdsgnmatshape,
		double *temppenmatshape, int *ntempcoeffshape, int *nppartempcoeffshape,
		double *temppenaltyshape, double *loccoeff, double *scalecoeff,
		double *shapecoeff, double *tempcoeffloc, double *tempcoeffscale,
		double *tempcoeffshape, double *dns){

  /* This function computes the log-likelihood for a spatial GEV model
     i.e. the GEV paramters are defined through a response surface
     that could be either is simple linear model or a p-spline. */

  int flag = usetempcov[0] + usetempcov[1] + usetempcov[2];
  double *trendlocs, *trendscales, *trendshapes;

  double *locs = (double *)R_alloc(*nSite, sizeof(double)),
    *scales = (double *)R_alloc(*nSite, sizeof(double)),
    *shapes = (double *)R_alloc(*nSite, sizeof(double));

  //Stage 1: Computing the GEV parameters using the design matrix
  *dns = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff,
		       shapecoeff, *nSite, *nloccoeff, *nscalecoeff, *nshapecoeff, locs,
		       scales, shapes);

  if (flag){
    trendlocs = (double *) R_alloc(*nObs, sizeof(double));
    trendscales = (double *) R_alloc(*nObs, sizeof(double));
    trendshapes = (double *) R_alloc(*nObs, sizeof(double));

    dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, tempcoeffloc,
		      tempcoeffscale, tempcoeffshape, *nSite, *nObs, usetempcov, *ntempcoeffloc,
		      *ntempcoeffscale, *ntempcoeffshape, trendlocs, trendscales, trendshapes);

    for (int i=0;i<*nSite;i++)
      for (int j=0;j<*nObs;j++)
	if (((scales[i] + trendscales[j]) <= 0) || ((shapes[i] + trendshapes[j]) <= -1)){
	  *dns = MINF;
	  return;
	}
  }

  else if (*dns != 0.0)
    return;

  //Stage 2: Compute the log-likelihood (assuming independence between
  //stations)
  if (flag){
    for (int i=0;i<*nSite;i++){
      for (int j=0;j<*nObs;j++){
	if (ISNA(data[j + i * *nObs]))
	  continue;

	if (fabs(shapes[i] + trendshapes[j]) <= 1e-6){
	  data[j + i * *nObs] = (data[j + i * *nObs] - locs[i] - trendlocs[j]) /
	    (scales[i] + trendscales[j]);

	  *dns += -log(scales[i] + trendscales[j]) - data[j + i * *nObs] -
	    exp(-data[j + i * *nObs]);
	}

      else{
	data[j + i * *nObs] = 1 + (shapes[i] + trendshapes[j]) *
	  (data[j + i * *nObs] - locs[i] - trendlocs[j]) / (scales[i] + trendscales[j]);

	  if (data[j + i * *nObs] <= 0){
	    *dns = MINF;
	    return;
	  }

	  *dns += -log(scales[i] + trendscales[j]) -
	    R_pow(data[j + i * *nObs], -1 / (shapes[i] + trendshapes[j])) -
	    (1 / (shapes[i] + trendshapes[j]) + 1) * log(data[j + i * *nObs]);
	}
      }
    }
  }

  else {
    for (int i=0;i<*nSite;i++){

      if (fabs(shapes[i]) <= 1e-6){
	for (int j=0;j<*nObs;j++){
	  if (ISNA(data[j + i * *nObs]))
	    continue;

	  data[j + i * *nObs] = (data[j + i * *nObs] - locs[i]) /
	    scales[i];

	  *dns += -log(scales[i]) - data[j + i * *nObs] -
	    exp(-data[j + i * *nObs]);
	}
      }

      else{
	for (int j=0;j<*nObs;j++){
	  if (ISNA(data[j + i * *nObs]))
	    continue;

	  data[j + i * *nObs] = 1 + shapes[i] * (data[j + i * *nObs] - locs[i]) /
	    scales[i];

	  if (data[j + i * *nObs] <= 0){
	    *dns = MINF;
	    return;
	  }

	  *dns += -log(scales[i]) - R_pow(data[j + i * *nObs], -1/shapes[i]) -
	    (1/shapes[i] + 1) * log(data[j + i * *nObs]);
	}
      }
    }
  }

  //Stage 3: Penalizing the llik if p-splines are used
  // 1- For the location parameter
  if (*locpenalty > 0)
    *dns -= penalization(locpenmat, loccoeff, *locpenalty,
			 *nloccoeff, *npparloc);

  // 2- For the scale parameter
  if (*scalepenalty > 0)
    *dns -= penalization(scalepenmat, scalecoeff, *scalepenalty,
			 *nscalecoeff, *npparscale);

  // 3- For the shape parameter
  if (*shapepenalty > 0)
    *dns -= penalization(shapepenmat, shapecoeff, *shapepenalty,
			 *nshapecoeff, *npparshape);

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

  return;
}


