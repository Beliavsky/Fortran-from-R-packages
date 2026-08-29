#include "header.h"

/* This file collects some pieces of code that are shared across
   different max-stable models. Indeed, the smith, geometric and
   brown-resnick models have the same bivariate distributions. */

void marginalPartSmith(int *start, int *nObs, int *nSite, double *data, double *frech,
		       double *mahalDist, double *locs, double *scales, double *shapes,
		       double *trendlocs, double *trendscales, double *trendshapes,
		       int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
		       int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
		       double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
		       double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
		       double *hess, double *grad){

  /* This function computes the gradient of the log-pairwise
     likelihood for which the bivariate CDF is similar to the one of
     the Smith model. Please note that only the marginal part is
     computed !!!! */

  const int nPairs = *nSite * (*nSite - 1) / 2;
  int currentPair = -1, i;

  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){
      currentPair++;

      int k;
      double imahal = 1 / mahalDist[currentPair], imahalSquare = imahal * imahal;

      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  int l, idx = *start;

	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    for (l=0;l<(*nloccoeff + *nscalecoeff + *nshapecoeff + *ntemploccoeff + *ntempscalecoeff + *ntempshapecoeff);l++)
	      hess[((idx + l) * *nObs + k) * nPairs + currentPair] = NA_REAL;

	    continue;
	  }

	  double ifrech1 = 1 / frech[k + i * *nObs], ifrech2 = 1 / frech[k + j * *nObs],
	    ifrech1Square = ifrech1 * ifrech1, ifrech2Square = ifrech2 * ifrech2,
	    c1 = log(frech[k + j * *nObs] * ifrech1) * imahal + 0.5 * mahalDist[currentPair],
	    c2 = mahalDist[currentPair] - c1,
	    dnormc1 = dnorm(c1, 0., 1., 0), pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0), pnormc2 = pnorm(c2, 0., 1., 1, 0),
	    dE, dz1loc, dz2loc, dz1scale, dz2scale, dz1shape, dz2shape;

	  double loc1 = locs[i] + trendlocs[k], loc2 = locs[j] + trendlocs[k],
	    scale1 = scales[i] + trendscales[k], scale2 = scales[j] + trendscales[k],
	    shape1 = shapes[i] + trendshapes[k], shape2 = shapes[j] + trendshapes[k];

	  //A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	  double B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	    dnormc2 * imahal * ifrech2Square,
	    C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	    dnormc1 * imahal * ifrech1Square,
	    D = c2 * dnormc1 * imahalSquare * ifrech1Square * ifrech2 +
	    c1 * dnormc2 * imahalSquare * ifrech1 * ifrech2Square,
	    iBCplusD = 1 / (B * C + D),
	    dAz1 = C,
	    dAz2 = B,
	    dBz1 = D,
	    dBz2 = (mahalDist[currentPair] + c1) * dnormc1 * imahalSquare *
	    ifrech1 * ifrech2Square - 2 * pnormc2 * ifrech2 * ifrech2Square -
	    (2 * mahalDist[currentPair] + c1) * dnormc2 * imahalSquare * ifrech2 *
	    ifrech2Square,
	    dCz1 = (mahalDist[currentPair] + c2) * dnormc2 * imahalSquare *
	    ifrech1Square * ifrech2 - 2 * pnormc1 * ifrech1Square * ifrech1 -
	    (2 * mahalDist[currentPair] + c2) * dnormc1 * imahalSquare * ifrech1 *
	    ifrech1Square,
	    dCz2 = D,
	    dDz1 = (1 - c2 * (mahalDist[currentPair] + c2)) * dnormc1 *
	    imahal * imahalSquare * ifrech1 * ifrech1Square * ifrech2 -
	    (1 + c1 * (mahalDist[currentPair] + c2)) * dnormc2 *
	    imahal * imahalSquare * ifrech1Square * ifrech2Square,
	    dDz2 = (1 - c1 * (mahalDist[currentPair] + c1)) * dnormc2 *
	    imahal * imahalSquare * ifrech1 * ifrech2 * ifrech2Square -
	    (1 + c2 * (mahalDist[currentPair] + c1)) * dnormc1 *
	    imahal * imahalSquare * ifrech1Square * ifrech2Square;

	  for (l=0;l<*nloccoeff;l++){
	    dE = (shape1 - 1) * locdsgnmat[i + *nSite * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      (shape2 - 1) * locdsgnmat[j + *nSite * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shape1) / scale1 * locdsgnmat[i + *nSite * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shape2) / scale2 * locdsgnmat[j + *nSite * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1loc + dAz2 * dz2loc) + ((dBz1 * dz1loc + dBz2 * dz2loc) * C +
						  B * (dCz1 * dz1loc + dCz2 * dz2loc) +
						  (dDz1 * dz1loc + dDz2 * dz2loc)) *
	       iBCplusD + dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];

	  }

	  idx += *nloccoeff;

	  for (l=0;l<*nscalecoeff;l++){
	    dE = scaledsgnmat[i + *nSite * l] * (loc1 - scale1 - data[k + i * *nObs]) /
	      (scale1 *  scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      scaledsgnmat[j + *nSite * l] * (loc2 - scale2 - data[k + j * *nObs]) /
	      (scale2 * scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shape1) * (data[k + i * *nObs] - loc1) /
	      (scale1 * scale1) * scaledsgnmat[i + *nSite * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shape2) * (data[k + j * *nObs] - loc2) /
	      (scale2 * scale2) * scaledsgnmat[j + *nSite * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1scale + dAz2 * dz2scale) + ((dBz1 * dz1scale + dBz2 * dz2scale) * C +
						      B * (dCz1 * dz1scale + dCz2 * dz2scale) +
						      (dDz1 * dz1scale + dDz2 * dz2scale)) *
	       iBCplusD + dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nscalecoeff;

	  for (l=0;l<*nshapecoeff;l++){
	    dE = -shapedsgnmat[i + *nSite * l] * log(frech[k + i * *nObs]) / shape1 + (1 / shape1 - 1) *
	      (data[k + i * *nObs] - loc1) * shapedsgnmat[i + *nSite * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) -
	      shapedsgnmat[j + *nSite * l] * log(frech[k + j * *nObs]) / shape2 + (1 / shape2 - 1) *
	      (data[k + j * *nObs] - loc2) * shapedsgnmat[j + *nSite * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1shape = frech[k + i * *nObs] * shapedsgnmat[i + *nSite * l] *
	      (-log(frech[k + i * *nObs]) / shape1 + (data[k + i * *nObs] - loc1) /
	       (shape1 * scale1 * R_pow(frech[k + i * *nObs], shape1)));
	    dz2shape = frech[k + j * *nObs] * shapedsgnmat[j + *nSite * l] *
	      (-log(frech[k + j * *nObs]) / shape2 + (data[k + j * *nObs] - loc2) /
	       (shape2 * scale2 * R_pow(frech[k + j * *nObs], shape2)));

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1shape + dAz2 * dz2shape) + ((dBz1 * dz1shape + dBz2 * dz2shape) * C +
						      B * (dCz1 * dz1shape + dCz2 * dz2shape) +
						      (dDz1 * dz1shape + dDz2 * dz2shape)) *
	       iBCplusD + dE);
	    grad[(idx + l) * *nObs + k] +=  hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nshapecoeff;

	  for (l=0;l<*ntemploccoeff;l++){
	    dE = (shape1 - 1) * tempdsgnmatloc[k + *nObs * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      (shape2 - 1) * tempdsgnmatloc[k + *nObs * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shape1) / scale1 * tempdsgnmatloc[k + *nObs * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shape2) / scale2 * tempdsgnmatloc[k + *nObs * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1loc + dAz2 * dz2loc) + ((dBz1 * dz1loc + dBz2 * dz2loc) * C +
						  B * (dCz1 * dz1loc + dCz2 * dz2loc) +
						  (dDz1 * dz1loc + dDz2 * dz2loc)) *
	       iBCplusD + dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *ntemploccoeff;

	  for (l=0;l<*ntempscalecoeff;l++){
	    dE = tempdsgnmatscale[k + *nObs * l] * (loc1 - scale1 - data[k + i * *nObs]) /
	      (scale1 *  scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      tempdsgnmatscale[k + *nObs * l] * (loc2 - scale2 - data[k + j * *nObs]) /
	      (scale2 * scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shape1) * (data[k + i * *nObs] - loc1) /
	      (scale1 * scale1) * tempdsgnmatscale[k + *nObs * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shape2) * (data[k + j * *nObs] - loc2) /
	      (scale2 * scale2) * tempdsgnmatscale[k + *nObs * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1scale + dAz2 * dz2scale) + ((dBz1 * dz1scale + dBz2 * dz2scale) * C +
						      B * (dCz1 * dz1scale + dCz2 * dz2scale) +
						      (dDz1 * dz1scale + dDz2 * dz2scale)) *
	       iBCplusD + dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *ntempscalecoeff;

	  for (l=0;l<*ntempshapecoeff;l++){
	    dE = -tempdsgnmatshape[k + *nObs * l] * log(frech[k + i * *nObs]) / shape1 +
	      (1 / shape1 - 1) * (data[k + i * *nObs] - loc1) * tempdsgnmatshape[k + *nObs * l] /
	      (scale1 * R_pow(frech[k + i * *nObs], shape1)) -
	      tempdsgnmatshape[k + *nObs * l] * log(frech[k + j * *nObs]) / shape2 +
	      (1 / shape2 - 1) * (data[k + j * *nObs] - loc2) * tempdsgnmatshape[k + *nObs * l] /
	      (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1shape = frech[k + i * *nObs] * tempdsgnmatshape[k + *nObs * l] *
	      (-log(frech[k + i * *nObs]) / shape1 + (data[k + i * *nObs] - loc1) /
	       (shape1 * scale1 * R_pow(frech[k + i * *nObs], shape1)));
	    dz2shape = frech[k + j * *nObs] * tempdsgnmatshape[k + *nObs * l] *
	      (-log(frech[k + j * *nObs]) / shape2 + (data[k + j * *nObs] - loc2) /
	       (shape2 * scale2 * R_pow(frech[k + j * *nObs], shape2)));

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1shape + dAz2 * dz2shape) + ((dBz1 * dz1shape + dBz2 * dz2shape) * C +
						      B * (dCz1 * dz1shape + dCz2 * dz2shape) +
						      (dDz1 * dz1shape + dDz2 * dz2shape)) *
	       iBCplusD + dE);
	    grad[(idx + l) * *nObs + k] +=  hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }

  return;

}

void marginalPartSchlat(int *start, int *nObs, int *nSite, double *data, double *frech,
			double *rho, double *locs, double *scales, double *shapes,
			double *trendlocs, double *trendscales, double *trendshapes,
			int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
			int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
			double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
			double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
			double *hess, double *grad){

  /* This function computes the gradient of the pairwise likelihood
     for the Schlather model. Note that only the marginal part is
     computed. */

  const int nPairs = *nSite * (*nSite - 1) / 2;
  int currentPair = -1, i;

  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){

      currentPair++;

      int k;
      double oneMinusRhoSquare = 1 - rho[currentPair] * rho[currentPair];

      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  int l, idx = *start;

	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    for (l=0;l<(*nloccoeff + *nscalecoeff + *nshapecoeff + *ntemploccoeff + *ntempscalecoeff + *ntempshapecoeff);l++)
	      hess[((idx + l) * *nObs + k) * nPairs + currentPair] = NA_REAL;

	    continue;
	  }


	  double frech1Square = frech[k + i * *nObs] * frech[k + i * *nObs],
	    frech2Square = frech[k + j * *nObs] * frech[k + j * *nObs];

	  double loc1 = locs[i] + trendlocs[k], loc2 = locs[j] + trendlocs[k],
	    scale1 = scales[i] + trendscales[k], scale2 = scales[j] + trendscales[k],
	    shape1 = shapes[i] + trendshapes[k], shape2 = shapes[j] + trendshapes[k],
	    dE, dz1loc, dz2loc, dz1scale, dz2scale, dz1shape, dz2shape;

	  double c1 = sqrt(frech1Square + frech2Square - 2 * frech[k + i * *nObs] *
			   frech[k + j * *nObs] * rho[currentPair]),
	    B = oneMinusRhoSquare / (2 * c1 * c1 * c1),
	    C = (frech[k + j * *nObs] + c1 - rho[currentPair] * frech[k + i * *nObs]) /
	    (2 * c1 * frech1Square),
	    D = (frech[k + i * *nObs] + c1 - rho[currentPair] * frech[k + j * *nObs]) /
	    (2 * c1 * frech2Square),
	    dAz1 = C,
	    dAz2 = D,
	    dBz1 = -3 * oneMinusRhoSquare * (frech[k + i * *nObs] - rho[currentPair] * frech[k + j * *nObs]) /
	    (2 * c1 * c1 * c1 * c1 * c1),
	    dBz2 = -3 * oneMinusRhoSquare * (frech[k + j * *nObs] - rho[currentPair] * frech[k + i * *nObs]) /
	    (2 * c1 * c1 * c1 * c1 * c1),
	    dCz1 = (2 * rho[currentPair] * frech1Square * frech[k + i * *nObs] + 6 * frech[k + i * *nObs] *
		    frech2Square * rho[currentPair] - 3 * frech1Square * frech[k + j * *nObs] *
		    (1 + rho[currentPair] * rho[currentPair]) - 2 * c1 * c1 * c1 - 2 * frech2Square *
		    frech[k + j * *nObs]) / (2 * c1 * c1 * c1 * frech1Square * frech[k + i * *nObs]),
	    dCz2 = B,
	    dDz1 = B,
	    dDz2 = (2 * rho[currentPair] * frech2Square * frech[k + j * *nObs] + 6 * frech[k + j * *nObs] *
		    frech1Square * rho[currentPair] - 3 * frech2Square * frech[k + i * *nObs] *
		    (1 + rho[currentPair] * rho[currentPair]) - 2 * c1 * c1 * c1 - 2 * frech1Square *
		    frech[k + i * *nObs]) / (2 * c1 * c1 * c1 * frech2Square * frech[k + j * *nObs]);

	  for (l=0;l<*nloccoeff;l++){
	    dE = (shape1 - 1) * locdsgnmat[i + *nSite * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      (shape2 - 1) * locdsgnmat[j + *nSite * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));
	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shape1) / scale1 * locdsgnmat[i + *nSite * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shape2) / scale2 * locdsgnmat[j + *nSite * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1loc + dAz2 * dz2loc) + ((dBz1 * dz1loc + dBz2 * dz2loc) +
						  (dCz1 * dz1loc + dCz2 * dz2loc) * D +
						  (dDz1 * dz1loc + dDz2 * dz2loc) * C) /
	       (B + C * D) + dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nloccoeff;

	  for (l=0;l<*nscalecoeff;l++){
	    dE = scaledsgnmat[i + *nSite * l] * (loc1 - scale1 - data[k + i * *nObs]) /
	      (scale1 * scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      scaledsgnmat[j + *nSite * l] * (loc2 - scale2 - data[k + j * *nObs]) /
	      (scale2 * scale2 * R_pow(frech[k + j * *nObs], shape2));
	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shape1) * (data[k + i * *nObs] - loc1) /
	      (scale1 * scale1) * scaledsgnmat[i + *nSite * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shape2) * (data[k + j * *nObs] - loc2) /
	      (scale2 * scale2) * scaledsgnmat[j + *nSite * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1scale + dAz2 * dz2scale) + ((dBz1 * dz1scale + dBz2 * dz2scale) +
						      (dCz1 * dz1scale + dCz2 * dz2scale) * D +
						      (dDz1 * dz1scale + dDz2 * dz2scale) * C) /
	       (B + C * D) + dE);
	    grad[(idx + l) * *nObs + k] +=  hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nscalecoeff;

	  for (l=0;l<*nshapecoeff;l++){
	    dE = -shapedsgnmat[i + *nSite * l] * log(frech[k + i * *nObs]) /
	      shape1 + (1/shape1 - 1) * (data[k + i * *nObs] - loc1) *
	      shapedsgnmat[i + *nSite * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) -
	      shapedsgnmat[j + *nSite * l] * log(frech[k + j * *nObs]) /
	      shape2 + (1/shape2 - 1) * (data[k + j * *nObs] - loc2) *
	      shapedsgnmat[j + *nSite * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));
	    dz1shape = frech[k + i * *nObs] * shapedsgnmat[i + *nSite * l] *
	      (-log(frech[k + i * *nObs]) / shape1 + (data[k + i * *nObs] - loc1) /
	       (shape1 * scale1 * R_pow(frech[k + i * *nObs], shape1)));
	    dz2shape = frech[k + j * *nObs] * shapedsgnmat[j + *nSite * l] *
	      (-log(frech[k + j * *nObs]) / shape2 + (data[k + j * *nObs] - loc2) /
	       (shape2 * scale2 * R_pow(frech[k + j * *nObs], shape2)));

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1shape + dAz2 * dz2shape) + ((dBz1 * dz1shape + dBz2 * dz2shape) +
						      (dCz1 * dz1shape + dCz2 * dz2shape) * D +
						      (dDz1 * dz1shape + dDz2 * dz2shape) * C) /
	       (B + C * D) + dE);
	    grad[(idx + l) * *nObs + k] +=  hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nshapecoeff;

	  for (l=0;l<*ntemploccoeff;l++){
	    dE = (shape1 - 1) * tempdsgnmatloc[k + *nObs * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      (shape2 - 1) * tempdsgnmatloc[k + *nObs * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));
	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shape1) / scale1 * tempdsgnmatloc[k + *nObs * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shape2) / scale2 * tempdsgnmatloc[k + *nObs * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1loc + dAz2 * dz2loc) + ((dBz1 * dz1loc + dBz2 * dz2loc) +
						  (dCz1 * dz1loc + dCz2 * dz2loc) * D +
						  (dDz1 * dz1loc + dDz2 * dz2loc) * C) /
	       (B + C * D) + dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *ntemploccoeff;

	  for (l=0;l<*ntempscalecoeff;l++){
	    dE = tempdsgnmatscale[k + *nObs * l] * (loc1 - scale1 - data[k + i * *nObs]) /
	      (scale1 * scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      tempdsgnmatscale[k + *nObs * l] * (loc2 - scale2 - data[k + j * *nObs]) /
	      (scale2 * scale2 * R_pow(frech[k + j * *nObs], shape2));
	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shape1) * (data[k + i * *nObs] - loc1) /
	      (scale1 * scale1) * tempdsgnmatscale[k + *nObs * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shape2) * (data[k + j * *nObs] - loc2) /
	      (scale2 * scale2) * tempdsgnmatscale[k + *nObs * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1scale + dAz2 * dz2scale) + ((dBz1 * dz1scale + dBz2 * dz2scale) +
						      (dCz1 * dz1scale + dCz2 * dz2scale) * D +
						      (dDz1 * dz1scale + dDz2 * dz2scale) * C) /
	       (B + C * D) + dE);
	    grad[(idx + l) * *nObs + k] +=  hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *ntempscalecoeff;

	  for (l=0;l<*ntempshapecoeff;l++){
	    dE = -tempdsgnmatshape[k + *nObs * l] * log(frech[k + i * *nObs]) /
	      shape1 + (1/shape1 - 1) * (data[k + i * *nObs] - loc1) *
	      tempdsgnmatshape[k + *nObs * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) -
	      tempdsgnmatshape[k + *nObs * l] * log(frech[k + j * *nObs]) /
	      shape2 + (1/shape2 - 1) * (data[k + j * *nObs] - loc2) *
	      tempdsgnmatshape[k + *nObs * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));
	    dz1shape = frech[k + i * *nObs] * tempdsgnmatshape[k + *nObs * l] *
	      (-log(frech[k + i * *nObs]) / shape1 + (data[k + i * *nObs] - loc1) /
	       (shape1 * scale1 * R_pow(frech[k + i * *nObs], shape1)));
	    dz2shape = frech[k + j * *nObs] * tempdsgnmatshape[k + *nObs * l] *
	      (-log(frech[k + j * *nObs]) / shape2 + (data[k + j * *nObs] - loc2) /
	       (shape2 * scale2 * R_pow(frech[k + j * *nObs], shape2)));

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1shape + dAz2 * dz2shape) + ((dBz1 * dz1shape + dBz2 * dz2shape) +
						      (dCz1 * dz1shape + dCz2 * dz2shape) * D +
						      (dDz1 * dz1shape + dDz2 * dz2shape) * C) /
	       (B + C * D) + dE);
	    grad[(idx + l) * *nObs + k] +=  hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }

  return;
}

void marginalPartiSchlat(int *start, int *nObs, int *nSite, double *data, double *frech,
			 double *alpha, double *rho, double *locs, double *scales, double *shapes,
			 double *trendlocs, double *trendscales, double *trendshapes,
			 int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
			 int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
			 double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
			 double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
			 double *hess, double *grad){

  /* This function computes the gradient of the pairwise likelihood
     for the independent Schlather model. Note that only the marginal
     part is computed. */

  const int nPairs = *nSite * (*nSite - 1) / 2;
  int currentPair = -1, i;

  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){

      currentPair++;
      int k;
      double oneMinusRhoSquare = 1 - rho[currentPair] * rho[currentPair];

      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  int l, idx = *start;

	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    for (l=0;l<(*nloccoeff + *nscalecoeff + *nshapecoeff + *ntemploccoeff + *ntempscalecoeff + *ntempshapecoeff);l++)
	      hess[((idx + l) * *nObs + k) * nPairs + currentPair] = NA_REAL;

	    continue;
	  }

	  double loc1 = locs[i] + trendlocs[k], loc2 = locs[j] + trendlocs[k],
	    scale1 = scales[i] + trendscales[k], scale2 = scales[j] + trendscales[k],
	    shape1 = shapes[i] + trendshapes[k], shape2 = shapes[j] + trendshapes[k],
	    dE, dz1loc, dz2loc, dz1scale, dz2scale, dz1shape, dz2shape;

	  double frech1Square = frech[k + i * *nObs] * frech[k + i * *nObs],
	    frech2Square = frech[k + j * *nObs] * frech[k + j * *nObs],
	    c1 = sqrt(frech1Square + frech2Square - 2 * frech[k + i * *nObs] * frech[k + j * *nObs] *
		      rho[currentPair]),
	    B = (1 - *alpha) * oneMinusRhoSquare / (2 * c1 * c1 * c1),
	    C = (*alpha - 1) * (rho[currentPair] * frech[k + i * *nObs] - c1 - frech[k + j * *nObs]) /
	    (2 * c1 * frech1Square) + *alpha / frech1Square,
	    D = (*alpha - 1) * (rho[currentPair] * frech[k + j * *nObs] - c1 - frech[k + i * *nObs]) /
	    (2 * c1 * frech2Square) + *alpha / frech2Square,
	    dAz1 = C,
	    dAz2 = D,
	    dBz1 = 3 * (*alpha - 1) * oneMinusRhoSquare * (frech[k + i * *nObs] - rho[currentPair] * frech[k + j * *nObs]) /
	    (2 * c1 * c1 * c1 * c1 * c1),
	    dBz2 = 3 * (*alpha - 1) * oneMinusRhoSquare * (frech[k + j * *nObs] - rho[currentPair] * frech[k + i * *nObs]) /
	    (2 * c1 * c1 * c1 * c1 * c1),
	    dCz1 = (2 * rho[currentPair] * frech1Square * frech[k + i * *nObs] + 6 * frech[k + i * *nObs] *
		    frech2Square * rho[currentPair] - 3 * frech1Square * frech[k + j * *nObs] *
		    (1 + rho[currentPair] * rho[currentPair]) - 2 * c1 * c1 * c1 -
		    2 * frech2Square * frech[k + j * *nObs]) / (2 * c1 * c1 * c1 * frech1Square * frech[k + i * *nObs]) *
	    (1 - *alpha) - 2 * *alpha / (frech1Square * frech[k + i * *nObs]),
	    dCz2 = B,
	    dDz1 = B,
	    dDz2 = (2 * rho[currentPair] * frech2Square * frech[k + j * *nObs] + 6 * frech[k + j * *nObs] *
		    frech1Square * rho[currentPair] - 3 * frech2Square * frech[k + i * *nObs] *
		    (1 + rho[currentPair] * rho[currentPair]) - 2 * c1 * c1 * c1 -
		    2 * frech1Square * frech[k + i * *nObs]) / (2 * c1 * c1 * c1 * frech2Square * frech[k + j * *nObs]) *
	    (1 - *alpha) - 2 * *alpha / (frech2Square * frech[k + j * *nObs]);

	  for (l=0;l<*nloccoeff;l++){
	    dE = (shape1 - 1) * locdsgnmat[i + *nSite * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      (shape2 - 1) * locdsgnmat[j + *nSite * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shape1) / scale1 * locdsgnmat[i + *nSite * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shape2) / scale2 * locdsgnmat[j + *nSite * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1loc + dAz2 * dz2loc) + ((dBz1 * dz1loc + dBz2 * dz2loc) +
						  (dCz1 * dz1loc + dCz2 * dz2loc) * D + C *
						  (dDz1 * dz1loc + dDz2 * dz2loc)) / (B + C * D) +
	       dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nloccoeff;

	  for (l=0;l<*nscalecoeff;l++){
	    dE = scaledsgnmat[i + *nSite * l] * (loc1 - scale1 - data[k + i * *nObs]) /
	      (scale1 * scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      scaledsgnmat[j + *nSite * l] * (loc2 - scale2 - data[k + j * *nObs]) /
	      (scale2 * scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shape1) * (data[k + i * *nObs] - loc1) /
	      (scale1 * scale1) * scaledsgnmat[i + *nSite * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shape2) * (data[k + j * *nObs] - loc2) /
	      (scale2 * scale2) * scaledsgnmat[j + *nSite * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1scale + dAz2 * dz2scale) + ((dBz1 * dz1scale + dBz2 * dz2scale) +
						      (dCz1 * dz1scale + dCz2 * dz2scale) * D +
						      C * (dDz1 * dz1scale + dDz2 * dz2scale)) /
	       (B + C * D) + dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nscalecoeff;

	  for (l=0;l<*nshapecoeff;l++){
	    dE = -shapedsgnmat[i + *nSite * l] * log(frech[k + i * *nObs]) / shape1 +
	      (1/shape1 - 1) * (data[k + i * *nObs] - loc1) * shapedsgnmat[i + *nSite * l] /
	      (scale1 * R_pow(frech[k + i * *nObs], shape1)) - shapedsgnmat[j + *nSite * l] *
	      log(frech[k + j * *nObs]) / shape2 + (1/shape2 - 1) * (data[k + j * *nObs] - loc2) *
	      shapedsgnmat[j + *nSite * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1shape = frech[k + i * *nObs] * shapedsgnmat[i + *nSite * l] *
	      (-log(frech[k + i * *nObs]) / shape1 + (data[k + i * *nObs] - loc1) /
	       (shape1 * scale1 * R_pow(frech[k + i * *nObs], shape1)));
	    dz2shape = frech[k + j * *nObs] * shapedsgnmat[j + *nSite * l] *
	      (-log(frech[k + j * *nObs]) / shape2 + (data[k + j * *nObs] - loc2) /
	       (shape2 * scale2 * R_pow(frech[k + j * *nObs], shape2)));

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1shape + dAz2 * dz2shape) + ((dBz1 * dz1shape + dBz2 * dz2shape) +
						      (dCz1 * dz1shape + dCz2 * dz2shape) * D +
						      (dDz1 * dz1shape + dDz2 * dz2shape) * C) /
	       (B + C * D) + dE);
	    grad[(idx + l) * *nObs + k] +=  hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nshapecoeff;

	  for (l=0;l<*ntemploccoeff;l++){
	    dE = (shape1 - 1) * tempdsgnmatloc[k + *nObs * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      (shape2 - 1) * locdsgnmat[j + *nSite * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shape1) / scale1 * tempdsgnmatloc[k + *nObs * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shape2) / scale2 * tempdsgnmatloc[k + *nObs * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1loc + dAz2 * dz2loc) + ((dBz1 * dz1loc + dBz2 * dz2loc) +
						  (dCz1 * dz1loc + dCz2 * dz2loc) * D + C *
						  (dDz1 * dz1loc + dDz2 * dz2loc)) / (B + C * D) +
	       dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *ntemploccoeff;

	  for (l=0;l<*ntempscalecoeff;l++){
	    dE = tempdsgnmatscale[k + *nObs * l] * (loc1 - scale1 - data[k + i * *nObs]) /
	      (scale1 * scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      tempdsgnmatscale[k + *nObs * l] * (loc2 - scale2 - data[k + j * *nObs]) /
	      (scale2 * scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shape1) * (data[k + i * *nObs] - loc1) /
	      (scale1 * scale1) * tempdsgnmatscale[k + *nObs * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shape2) * (data[k + j * *nObs] - loc2) /
	      (scale2 * scale2) * tempdsgnmatscale[j + *nObs * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1scale + dAz2 * dz2scale) + ((dBz1 * dz1scale + dBz2 * dz2scale) +
						      (dCz1 * dz1scale + dCz2 * dz2scale) * D +
						      C * (dDz1 * dz1scale + dDz2 * dz2scale)) /
	       (B + C * D) + dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *ntempscalecoeff;

	  for (l=0;l<*ntempshapecoeff;l++){
	    dE = -tempdsgnmatshape[k + *nObs * l] * log(frech[k + i * *nObs]) / shape1 +
	      (1/shape1 - 1) * (data[k + i * *nObs] - loc1) * tempdsgnmatshape[k + *nObs * l] /
	      (scale1 * R_pow(frech[k + i * *nObs], shape1)) - tempdsgnmatshape[k + *nObs * l] *
	      log(frech[k + j * *nObs]) / shape2 + (1/shape2 - 1) * (data[k + j * *nObs] - loc2) *
	      tempdsgnmatshape[k + *nObs * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1shape = frech[k + i * *nObs] * tempdsgnmatshape[k + *nObs * l] *
	      (-log(frech[k + i * *nObs]) / shape1 + (data[k + i * *nObs] - loc1) /
	       (shape1 * scale1 * R_pow(frech[k + i * *nObs], shape1)));
	    dz2shape = frech[k + j * *nObs] * tempdsgnmatshape[k + *nObs * l] *
	      (-log(frech[k + j * *nObs]) / shape2 + (data[k + j * *nObs] - loc2) /
	       (shape2 * scale2 * R_pow(frech[k + j * *nObs], shape2)));

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1shape + dAz2 * dz2shape) + ((dBz1 * dz1shape + dBz2 * dz2shape) +
						      (dCz1 * dz1shape + dCz2 * dz2shape) * D +
						      (dDz1 * dz1shape + dDz2 * dz2shape) * C) /
	       (B + C * D) + dE);
	    grad[(idx + l) * *nObs + k] +=  hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }

  return;
}

void marginalPartExtremalt(int *start, int *nObs, int *nSite, double *data, double *frech,
			   double *df, double *rho, double *locs, double *scales, double *shapes,
			   double *trendlocs, double *trendscales, double *trendshapes,
			   int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
			   int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
			   double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
			   double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
			   double *hess, double *grad){

  /* This function computes the gradient of the pairwise likelihood
     for the Extremal t model. Note that only the marginal part is
     computed. */

  const int nPairs = *nSite * (*nSite - 1) / 2;
  const double idf = 1 / *df, dfPlus1 = *df + 1;
  int currentPair = -1, i;

  for (i=0;i<(*nSite-1);i++){
    int j;
    for (j=i+1;j<*nSite;j++){

      currentPair++;
      int k;
      double a = sqrt(dfPlus1 / (1 - rho[currentPair] * rho[currentPair]));

      if (weights[currentPair] != 0){
	for (k=*nObs;k--;){
	  int l, idx = *start;

	  if (ISNA(frech[k + i * *nObs]) || ISNA(frech[k + j * *nObs])){
	    for (l=0;l<(*nloccoeff + *nscalecoeff + *nshapecoeff + *ntemploccoeff + *ntempscalecoeff + *ntempshapecoeff);l++)
	      hess[((idx + l) * *nObs + k) * nPairs + currentPair] = NA_REAL;

	    continue;
	  }

	  double loc1 = locs[i] + trendlocs[k], loc2 = locs[j] + trendlocs[k],
	    scale1 = scales[i] + trendscales[k], scale2 = scales[j] + trendscales[k],
	    shape1 = shapes[i] + trendshapes[k], shape2 = shapes[j] + trendshapes[k],
	    dE, dz1loc, dz2loc, dz1scale, dz2scale, dz1shape, dz2shape;

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
	    dertc1 = -(1 + 1 / dfPlus1) * c1 / (1  + c1 * c1 / dfPlus1) * tc1,
	    dertc2 = -(1 + 1 / dfPlus1) * c2 / (1  + c2 * c2 / dfPlus1) * tc2,
	    der2tc1 = dertc1 * dertc1 / tc1 + dertc1 / c1 + 2 * dertc1 * dertc1 /
	    (tc1 * (dfPlus1 + 1)),
	    der2tc2 = dertc2 * dertc2 / tc2 + dertc2 / c2 + 2 * dertc2 * dertc2 /
	    (tc2 * (dfPlus1 + 1)),
	    //partial derivatives of frechi_j w.r.t. z1 and z2
	    dfrech12_z1 = frech1_2 * idf * ifrech1,
	    dfrech12_z2 = -frech1_2 * idf * ifrech2,
	    dfrech21_z1 = -frech2_1 * idf * ifrech1,
	    dfrech21_z2 = frech2_1 * idf * ifrech2,
	    //partial derivatives of c1, c2 w.r.t. z1 and z2
	    dc1_z1 = a * dfrech21_z1,
	    dc1_z2 = a * dfrech21_z2,
	    dc2_z1 = a * dfrech12_z1,
	    dc2_z2 = a * dfrech12_z2,
	    //A = -Tc1 * ifrech1 - Tc2 * ifrech2,
	    B = ifrech1 * ifrech1 * Tc1 + ifrech1 * ifrech1 * idf * frech2_1 * tc1 * a -
	    ifrech1 * ifrech2 * idf * frech1_2 * tc2 * a,
	    C = ifrech2 * ifrech2 * Tc2 + ifrech2 * ifrech2 * idf * frech1_2 * tc2 * a -
	    ifrech1 * ifrech2 * idf * frech2_1 * tc1 * a,
	    D = ifrech1 * ifrech1 * ifrech2 * idf * idf * dfPlus1 * frech2_1 * a * tc1 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * dfPlus1 * frech1_2 * a * tc2 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a * dertc1 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a * dertc2,
	    iBCplusD = 1 / (B * C + D),
	    dAz1 = B,
	    dAz2 = C,
	    /* B = ifrech1 * ifrech1 * Tc1 + ifrech1 * ifrech1 * idf * frech2_1 * tc1 * a -
	       ifrech1 * ifrech2 * idf * frech1_2 * tc2 * a, */
	    dBz1 = -2 * ifrech1 * ifrech1 * ifrech1 * Tc1 + ifrech1 * ifrech1 * dc1_z1 * tc1 -
	    2 * ifrech1 * ifrech1 * ifrech1 * idf * frech2_1 * tc1 * a +
	    ifrech1 * ifrech1 * idf * dfrech21_z1 * tc1 * a +
	    ifrech1 * ifrech1 * idf * frech2_1 * dc1_z1 * dertc1 * a +
	    ifrech1 * ifrech1 * ifrech2 * idf * frech1_2 * tc2 * a -
	    ifrech1 * ifrech2 * idf * dfrech12_z1 * tc2 * a -
	    ifrech1 * ifrech2 * idf * frech1_2 * dc2_z1 * dertc2 * a,
	    dBz2 = D,
	    dCz1 = D,
	    /* C = ifrech2 * ifrech2 * Tc2 + ifrech2 * ifrech2 * idf * frech1_2 * tc2 * a -
	       ifrech1 * ifrech2 * idf * frech2_1 * tc1 * a, */
	    dCz2 = -2 * ifrech2 * ifrech2 * ifrech2 * Tc2 + ifrech2 * ifrech2 * dc2_z2 * tc2 -
	    2 * ifrech2 * ifrech2 * ifrech2 * idf * frech1_2 * tc2 * a +
	    ifrech2 * ifrech2 * idf * dfrech12_z2 * tc2 * a +
	    ifrech2 * ifrech2 * idf * frech1_2 * dc2_z2 * dertc2 * a +
	    ifrech1 * ifrech2 * ifrech2 * idf * frech2_1 * tc1 * a -
	    ifrech1 * ifrech2 * idf * dfrech21_z2 * tc1 * a -
	    ifrech1 * ifrech2 * idf * frech2_1 * dc1_z2 * dertc1 * a,
	    /* D = ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * dfPlus1 * a * tc1 +
	       ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * dfPlus1 * a * tc2 +
	       ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a * dertc1 +
	       ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a * dertc2, */
	    dDz1 = -2 * ifrech1 * ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * dfPlus1 *
	    a * tc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * dfPlus1 * dfrech21_z1 * a * tc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * dfPlus1 * frech2_1 * a * dc1_z1 * dertc1 -
	    ifrech1 * ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * dfPlus1 * a * tc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * dfPlus1 * dfrech12_z1 * a * tc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * dfPlus1 * frech1_2 * a * dc2_z1 * dertc2 -
	    2 * ifrech1 * ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a *
	    dertc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * 2 * dfrech21_z1 * frech2_1 * a * a *
	    dertc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a *
	    dc1_z1 * der2tc1 -
	    ifrech1 * ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a *
	    dertc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * 2 * dfrech12_z1 * frech1_2 * a * a *
	    dertc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a *
	    dc2_z1 * der2tc2,
	    /*D = ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * dfPlus1 * a * tc1 +
	      ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * dfPlus1 * a * tc2 +
	      ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a * dertc1 +
	      ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a * dertc2,*/
	    dDz2 = -ifrech1 * ifrech1 * ifrech2 * ifrech2 * idf * idf * frech2_1 * dfPlus1 *
	    a * tc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * dfrech21_z2 * dfPlus1 * a * tc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * dfPlus1 * a * dc1_z2 * dertc1 -
	    ifrech1 * 2 * ifrech2 * ifrech2 * ifrech2 * idf * idf * frech1_2 * dfPlus1 * a * tc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * dfrech12_z2 * dfPlus1 * a * tc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * dfPlus1 * a * dc2_z2 * dertc2 -
	    ifrech1 * ifrech1 * ifrech2 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a *
	    dertc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * 2 * dfrech21_z2 * frech2_1 * a * a *
	    dertc1 +
	    ifrech1 * ifrech1 * ifrech2 * idf * idf * frech2_1 * frech2_1 * a * a * dc1_z2 *
	    der2tc1 -
	    ifrech1 * 2 * ifrech2 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a *
	    dertc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * 2 * dfrech12_z2 * frech1_2 * a * a *
	    dertc2 +
	    ifrech1 * ifrech2 * ifrech2 * idf * idf * frech1_2 * frech1_2 * a * a * dc2_z2 *
	    der2tc2;

	  for (l=0;l<*nloccoeff;l++){
	    dE = (shape1 - 1) * locdsgnmat[i + *nSite * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      (shape2 - 1) * locdsgnmat[j + *nSite * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shape1) / scale1 * locdsgnmat[i + *nSite * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shape2) / scale2 * locdsgnmat[j + *nSite * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1loc + dAz2 * dz2loc) + ((dDz1 * dz1loc + dDz2 * dz2loc) +
						  (dBz1 * dz1loc + dBz2 * dz2loc) * C + B *
						  (dCz1 * dz1loc + dCz2 * dz2loc)) * iBCplusD +
	       dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nloccoeff;

	  for (l=0;l<*nscalecoeff;l++){
	    dE = scaledsgnmat[i + *nSite * l] * (loc1 - scale1 - data[k + i * *nObs]) /
	      (scale1 * scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      scaledsgnmat[j + *nSite * l] * (loc2 - scale2 - data[k + j * *nObs]) /
	      (scale2 * scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shape1) * (data[k + i * *nObs] - loc1) /
	      (scale1 * scale1) * scaledsgnmat[i + *nSite * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shape2) * (data[k + j * *nObs] - loc2) /
	      (scale2 * scale2) * scaledsgnmat[j + *nSite * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1scale + dAz2 * dz2scale) + ((dDz1 * dz1scale + dDz2 * dz2scale) +
						      (dBz1 * dz1scale + dBz2 * dz2scale) * C +
						      B * (dCz1 * dz1scale + dCz2 * dz2scale)) *
	       iBCplusD + dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nscalecoeff;

	  for (l=0;l<*nshapecoeff;l++){
	    dE = -shapedsgnmat[i + *nSite * l] * log(frech[k + i * *nObs]) / shape1 +
	      (1/shape1 - 1) * (data[k + i * *nObs] - loc1) * shapedsgnmat[i + *nSite * l] /
	      (scale1 * R_pow(frech[k + i * *nObs], shape1)) - shapedsgnmat[j + *nSite * l] *
	      log(frech[k + j * *nObs]) / shape2 + (1/shape2 - 1) * (data[k + j * *nObs] - loc2) *
	      shapedsgnmat[j + *nSite * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1shape = frech[k + i * *nObs] * shapedsgnmat[i + *nSite * l] *
	      (-log(frech[k + i * *nObs]) / shape1 + (data[k + i * *nObs] - loc1) /
	       (shape1 * scale1 * R_pow(frech[k + i * *nObs], shape1)));
	    dz2shape = frech[k + j * *nObs] * shapedsgnmat[j + *nSite * l] *
	      (-log(frech[k + j * *nObs]) / shape2 + (data[k + j * *nObs] - loc2) /
	       (shape2 * scale2 * R_pow(frech[k + j * *nObs], shape2)));

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1shape + dAz2 * dz2shape) + ((dDz1 * dz1shape + dDz2 * dz2shape) +
						      (dBz1 * dz1shape + dBz2 * dz2shape) * C +
						      (dCz1 * dz1shape + dCz2 * dz2shape) * B) *
	       iBCplusD + dE);
	    grad[(idx + l) * *nObs + k] +=  hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *nshapecoeff;

	  for (l=0;l<*ntemploccoeff;l++){
	    dE = (shape1 - 1) * tempdsgnmatloc[k + *nObs * l] / (scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      (shape2 - 1) * locdsgnmat[j + *nSite * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shape1) / scale1 * tempdsgnmatloc[k + *nObs * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shape2) / scale2 * tempdsgnmatloc[k + *nObs * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1loc + dAz2 * dz2loc) + ((dDz1 * dz1loc + dDz2 * dz2loc) +
						  (dBz1 * dz1loc + dBz2 * dz2loc) * C + B *
						  (dCz1 * dz1loc + dCz2 * dz2loc)) * iBCplusD +
	       dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *ntemploccoeff;

	  for (l=0;l<*ntempscalecoeff;l++){
	    dE = tempdsgnmatscale[k + *nObs * l] * (loc1 - scale1 - data[k + i * *nObs]) /
	      (scale1 * scale1 * R_pow(frech[k + i * *nObs], shape1)) +
	      tempdsgnmatscale[k + *nObs * l] * (loc2 - scale2 - data[k + j * *nObs]) /
	      (scale2 * scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shape1) * (data[k + i * *nObs] - loc1) /
	      (scale1 * scale1) * tempdsgnmatscale[k + *nObs * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shape2) * (data[k + j * *nObs] - loc2) /
	      (scale2 * scale2) * tempdsgnmatscale[j + *nObs * l];

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1scale + dAz2 * dz2scale) + ((dDz1 * dz1scale + dDz2 * dz2scale) +
						      (dBz1 * dz1scale + dBz2 * dz2scale) * C +
						      B * (dCz1 * dz1scale + dCz2 * dz2scale)) *
	       iBCplusD + dE);
	    grad[(idx + l) * *nObs + k] += hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }

	  idx += *ntempscalecoeff;

	  for (l=0;l<*ntempshapecoeff;l++){
	    dE = -tempdsgnmatshape[k + *nObs * l] * log(frech[k + i * *nObs]) / shape1 +
	      (1/shape1 - 1) * (data[k + i * *nObs] - loc1) * tempdsgnmatshape[k + *nObs * l] /
	      (scale1 * R_pow(frech[k + i * *nObs], shape1)) - tempdsgnmatshape[k + *nObs * l] *
	      log(frech[k + j * *nObs]) / shape2 + (1/shape2 - 1) * (data[k + j * *nObs] - loc2) *
	      tempdsgnmatshape[k + *nObs * l] / (scale2 * R_pow(frech[k + j * *nObs], shape2));

	    dz1shape = frech[k + i * *nObs] * tempdsgnmatshape[k + *nObs * l] *
	      (-log(frech[k + i * *nObs]) / shape1 + (data[k + i * *nObs] - loc1) /
	       (shape1 * scale1 * R_pow(frech[k + i * *nObs], shape1)));
	    dz2shape = frech[k + j * *nObs] * tempdsgnmatshape[k + *nObs * l] *
	      (-log(frech[k + j * *nObs]) / shape2 + (data[k + j * *nObs] - loc2) /
	       (shape2 * scale2 * R_pow(frech[k + j * *nObs], shape2)));

	    hess[((idx + l) * *nObs + k) * nPairs + currentPair] = weights[currentPair] *
	      ((dAz1 * dz1shape + dAz2 * dz2shape) + ((dDz1 * dz1shape + dDz2 * dz2shape) +
						      (dBz1 * dz1shape + dBz2 * dz2shape) * C +
						      (dCz1 * dz1shape + dCz2 * dz2shape) * B) *
	       iBCplusD + dE);
	    grad[(idx + l) * *nObs + k] +=  hess[((idx + l) * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }

  return;
}

