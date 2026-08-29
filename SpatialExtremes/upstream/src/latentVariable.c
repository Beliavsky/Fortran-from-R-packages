#include "header.h"

void latentgev(int *n, double *data, int *nSite, int *nObs, int *covmod,
	       int *dim, double *distMat, double *dsgnMat, int *nBeta,
	       double *beta, double *sills, double *ranges, double *smooths,
	       double *gevParams, double *hyperSill, double *hyperRange,
	       double *hyperSmooth, double *hyperBetaMean,
	       double *hyperBetaIcov, double *propGev, double *propRanges,
	       double *propSmooths, int *useLogLink, double *mcLoc, double *mcScale,
	       double *mcShape, double *accRates, double *extRates, int *thin,
	       int *burnin){


  int iter = 0, iterThin = 0, idxSite, idxSite2, idxMarge, idxBeta, info = 0,
    oneInt = 1, nSite2 = *nSite * *nSite,
    nPairs = *nSite * (*nSite + 1) / 2,
    *cumBeta = (int *) R_alloc(4, sizeof(int)),
    *cumBeta2 = (int *) R_alloc(3, sizeof(int)),
    *nBeta2 = (int *) R_alloc(3, sizeof(int)),
    lagLoc = nBeta[0] + 3 + *nSite, lagScale = nBeta[1] + 3 + *nSite,
    lagShape = nBeta[2] + 3 + *nSite;

  cumBeta[0] = 0;
  cumBeta[1] = nBeta[0];
  cumBeta[2] = nBeta[0] + nBeta[1];
  cumBeta[3] = cumBeta[2] + nBeta[2];
  cumBeta2[0] = 0;
  cumBeta2[1] = nBeta[0] * nBeta[0];
  cumBeta2[2] = nBeta[0] * nBeta[0] + nBeta[1] * nBeta[1];
  nBeta2[0] = nBeta[0] * nBeta[0];
  nBeta2[1] = nBeta[1] * nBeta[1];
  nBeta2[2] = nBeta[2] * nBeta[2];

  double one = 1.0, zero = 0.0, flag = 0.0, logDetProp,
    *logDet = (double *) R_alloc(3, sizeof(double)),
    *covMatChol = (double *) R_alloc(3 * nSite2, sizeof(double)),
    *GPmean = (double *) R_alloc(3 * *nSite, sizeof(double)),
    *resTop = (double *) R_alloc(*nSite, sizeof(double)),
    *resBottom = (double *) R_alloc(*nSite, sizeof(double)),
    *covariances = (double *) R_alloc(nPairs, sizeof(double)),
    *proposalGEV = (double *) R_alloc(3, sizeof(double)),
    *covMatPropChol = (double *) R_alloc(nSite2, sizeof(double));

  for (int i=3;i--;)
    logDet[i] = 0;

  for (int i=(3 * nSite2);i--;)
    covMatChol[i] = 0;

  for (int i=(3 * *nSite);i--;)
    GPmean[i] = 0;

  for (int i=nSite2;i--;)
    covMatPropChol[i] = 0;

  /*----------------------------------------------------*/
  /*                                                    */
  /*           Compute some initial objects             */
  /*                                                    */
  /*----------------------------------------------------*/

  // a. The inverse of the covariance matrices
  for (idxMarge=0;idxMarge<3;idxMarge++){

    switch(covmod[idxMarge]){
    case 1:
      flag = whittleMatern(distMat, nPairs, zero, sills[idxMarge], ranges[idxMarge],
			   smooths[idxMarge], covariances);
      break;
    case 2:
      flag = cauchy(distMat, nPairs, zero, sills[idxMarge], ranges[idxMarge],
		    smooths[idxMarge], covariances);
      break;
    case 3:
      flag = powerExp(distMat, nPairs, zero, sills[idxMarge], ranges[idxMarge],
		      smooths[idxMarge], covariances);
      break;
    case 4:
      flag = bessel(distMat, nPairs, *dim, zero, sills[idxMarge], ranges[idxMarge],
		    smooths[idxMarge], covariances);
      break;
    }

    if (flag != 0)
      error("The starting values (covariance parameter) are ill-defined. Please check\n");

    /* We need to fill in the upper triangular part of covMatChol with
       covariances */
    {
      int current=-1;
      for (idxSite=0;idxSite<*nSite;idxSite++)
	for (idxSite2=idxSite;idxSite2<*nSite;idxSite2++){
	  current++;
	  covMatChol[idxSite + idxSite2 * *nSite + idxMarge * nSite2] = covariances[current];
	}
    }

    // Finally compute its Cholesky decomposition
    F77_CALL(dpotrf)("U", nSite, covMatChol + idxMarge * nSite2, nSite, &info FCONE);

    if (info != 0)
      error("Impossible to get the Cholesky decomp. from the starting values\n");

    /* Compute the log of the determinant of the proposal
       cov. mat. using the sum of the square of the diagonal elements of
       the Cholesky decomposition */
    for (idxSite2=0;idxSite2<*nSite;idxSite2++)
      logDet[idxMarge] += log(covMatChol[idxSite2 * (*nSite + 1) + idxMarge *
					 nSite2]);

    logDet[idxMarge] *= 2;
  }

  // b. The mean of the Gaussian processes
  for (idxMarge=0;idxMarge<3;idxMarge++)
    for (idxSite=0;idxSite<*nSite;idxSite++)
      for (idxBeta=0;idxBeta<nBeta[idxMarge];idxBeta++)
	GPmean[idxSite + idxMarge * *nSite] +=
	  dsgnMat[idxBeta * *nSite + idxSite + cumBeta[idxMarge] * *nSite] *
	  beta[cumBeta[idxMarge] + idxBeta];

  // c. Some constant related to the conjugate distributions
  double *conjMeanCst = (double *)R_alloc(cumBeta[3], sizeof(double));
  for(int i=cumBeta[3];i--;)
    conjMeanCst[i]=0;

  for (idxMarge=0;idxMarge<3;idxMarge++)
    F77_CALL(dsymv)("U", nBeta + idxMarge, &one, hyperBetaIcov +
		    cumBeta2[idxMarge], nBeta + idxMarge, hyperBetaMean +
		    cumBeta[idxMarge], &oneInt, &zero, conjMeanCst + cumBeta[idxMarge],
		    &oneInt FCONE);

  /*----------------------------------------------------*/
  /*                                                    */
  /*               Starting the MCMC algo               */
  /*                                                    */
  /*----------------------------------------------------*/

  GetRNGstate();
  while (iterThin<*n){
    void R_CheckUserInterrupt(void);

    /*----------------------------------------------------*/
    /*                                                    */
    /*           Updating the GEV parameters              */
    /*                                                    */
    /*----------------------------------------------------*/

    for (idxSite=0;idxSite<*nSite;idxSite++){
      for (idxMarge=0;idxMarge<3;idxMarge++){
	double logpropRatio = 0;
	proposalGEV[0] = gevParams[idxSite];
	proposalGEV[1] = gevParams[*nSite + idxSite];
	proposalGEV[2] = gevParams[2 * *nSite + idxSite];

	if (idxMarge==1){
	  proposalGEV[1] = rlnorm(log(gevParams[*nSite + idxSite]), propGev[1]);
	  logpropRatio = log(proposalGEV[1] / gevParams[*nSite + idxSite]);
	}

	else
	  proposalGEV[idxMarge] = rnorm(gevParams[idxMarge * *nSite + idxSite], propGev[idxMarge]);

	double topGEV = 0, bottomGEV = 0;
	gevlik(data + idxSite * *nObs, nObs, proposalGEV, proposalGEV + 1,
	       proposalGEV + 2, &topGEV);

	if (topGEV == -1e6){
	  extRates[idxMarge]++;
	  continue;
	}

	gevlik(data + idxSite * *nObs, nObs, gevParams + idxSite, gevParams +
	       *nSite + idxSite, gevParams + 2 * *nSite + idxSite, &bottomGEV);

	double topGP = 0, bottomGP = 0;
	if ((idxMarge==1) && (*useLogLink)){
	  for (idxSite2=0;idxSite2<*nSite;idxSite2++)
	    resBottom[idxSite2] = log(gevParams[idxSite2 + idxMarge * *nSite]) -
	      GPmean[idxSite2 + idxMarge * *nSite];
	  
	  memcpy(resTop, resBottom, *nSite * sizeof(double));
	  resTop[idxSite] = log(proposalGEV[idxMarge]) - GPmean[idxSite + idxMarge * *nSite];
	} else {
	  for (idxSite2=0;idxSite2<*nSite;idxSite2++)
	    resBottom[idxSite2] = gevParams[idxSite2 + idxMarge * *nSite] -
	      GPmean[idxSite2 + idxMarge * *nSite];

	  memcpy(resTop, resBottom, *nSite * sizeof(double));
	  resTop[idxSite] = proposalGEV[idxMarge] - GPmean[idxSite + idxMarge *
							   *nSite];
	}

	F77_CALL(dtrsm)("L", "U", "T", "N", nSite, &oneInt, &one, covMatChol +
			idxMarge * nSite2, nSite, resTop, nSite
			FCONE FCONE FCONE FCONE);
	F77_CALL(dtrsm)("L", "U", "T", "N", nSite, &oneInt, &one, covMatChol +
			idxMarge * nSite2, nSite, resBottom, nSite
			FCONE FCONE FCONE FCONE);

	for (idxSite2=0;idxSite2<*nSite;idxSite2++){
	  topGP += resTop[idxSite2] * resTop[idxSite2];
	  bottomGP += resBottom[idxSite2] * resBottom[idxSite2];
	}

	topGP *= -0.5;
	bottomGP *= -0.5;

	if ((idxMarge==1) && (*useLogLink)){
	  topGP -= log(proposalGEV[idxMarge]);
	  bottomGP -= log(gevParams[*nSite + idxSite]);
	}

	if (unif_rand() < exp(topGEV - bottomGEV + topGP - bottomGP +
			      logpropRatio)){
	  gevParams[idxSite + idxMarge * *nSite] = proposalGEV[idxMarge];
	  accRates[idxMarge]++;
	}
      }
    }

    /*----------------------------------------------------*/
    /*                                                    */
    /*        Updating the regression parameters          */
    /*                (conjugate prior)                   */
    /*                                                    */
    /*----------------------------------------------------*/

    for (idxMarge=0;idxMarge<3;idxMarge++){

      /* conjCovMat is the covariance matrix for the conjugate
	 distribution i.e. MVN

	 conjCovMatChol is its Cholesky decomposition */
      double *dummy = malloc(*nSite * nBeta[idxMarge] * sizeof(double)),
	*conjCovMat = malloc(nBeta2[idxMarge] * sizeof(double)),
	*conjCovMatChol = malloc(nBeta2[idxMarge] * sizeof(double));

      memcpy(conjCovMat, hyperBetaIcov + cumBeta2[idxMarge],
	     nBeta2[idxMarge] * sizeof(double));
      memcpy(dummy, dsgnMat + *nSite * cumBeta[idxMarge],
	     *nSite * nBeta[idxMarge] * sizeof(double));

      // Compute dummy = covMatChol^(-T) %*% dsgnMat
      F77_CALL(dtrsm)("L", "U", "T", "N", nSite, nBeta + idxMarge, &one,
		      covMatChol + idxMarge * nSite2, nSite, dummy, nSite
		      FCONE FCONE FCONE FCONE);

      /* Compute conjCovMat = dummy^T %*% dummy + conjCovMat

	 WARNING: Only the upper diagonal elements will be stored */
      F77_CALL(dsyrk)("U", "T", nBeta + idxMarge, nSite, &one, dummy, nSite,
		      &one, conjCovMat, nBeta + idxMarge FCONE FCONE);

      /* Rmk: The "real" conjugate cov. matrix is the inverse of
	 conjCovMat but it is not necessary to compute it */

      //Compute its Cholesky decomposition
      memcpy(conjCovMatChol, conjCovMat, nBeta2[idxMarge] * sizeof(double));
      F77_CALL(dpotrf)("U", nBeta + idxMarge, conjCovMatChol, nBeta + idxMarge,
		       &info FCONE);

      // Compute dummy2 = covMatChol^(-T) %*% (locs or scales or shapes)
      double *dummy2 = malloc(*nSite * sizeof(double));
      memcpy(dummy2, gevParams + idxMarge * *nSite, *nSite * sizeof(double));

      if ((idxMarge==1) && *useLogLink){
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  dummy2[idxSite] = log(dummy2[idxSite]);
      }
      
      F77_CALL(dtrsm)("L", "U", "T", "N", nSite, &oneInt, &one, covMatChol +
		      idxMarge * nSite2, nSite, dummy2, nSite
		      FCONE FCONE FCONE FCONE);

      // conjMean is the mean for the conjugate distribution i.e. MVN
      // Set conjMean = conjMeanCst := hyperBetaIcov %*% hyperBetaMean
      double *conjMean = malloc(nBeta[idxMarge] * sizeof(double));
      memcpy(conjMean, conjMeanCst + cumBeta[idxMarge],
	     nBeta[idxMarge] * sizeof(double));

      // Compute conjMean = conjMean + dummy^T %*% dummy2 (dummy2 is a vector)
      F77_CALL(dgemv)("T", nSite, nBeta + idxMarge, &one, dummy, nSite, dummy2,
		      &oneInt, &one, conjMean, &oneInt FCONE);

      // Compute conjMean = conjCovMat^(-1) %*% conjMean
      F77_CALL(dposv)("U", nBeta + idxMarge, &oneInt, conjCovMat, nBeta +
		      idxMarge, conjMean, nBeta + idxMarge, &info FCONE);

      /* The new state is a realisation from the MVN(conjMean,
	 conjCovMat) so we simulate it from the Cholesky
	 decomposition */

      double *stdNormal = malloc(nBeta[idxMarge] * sizeof(double));
      for (idxBeta=0;idxBeta<nBeta[idxMarge];idxBeta++)
	stdNormal[idxBeta] = norm_rand();

      /* Rmk: Recall that conjCovMat is the precision matrix and *NOT*
	 the covariance matrix. Instead of using the Cholesky
	 decomposition of the conjugate covariance matrix (that we
	 still haven't computed), we use the inverse of the Cholesky
	 decomposition. This is different from the standard simulation
	 technique but completely equivalent since

	      iSigma = iSigma_*^T %*% iSigma_*
	 <==> Sigma := iSigma^(-1) = iSigma_*^(-1) %*% iSigma_*^(-T),

	 where iSigma_* is the Cholesky decomposition of iSigma.

	 Therefore we can use iSigma_*^(-1) for the simulation. */
      F77_CALL(dtrsm)("L", "U", "N", "N", nBeta + idxMarge, &oneInt,
		      &one, conjCovMatChol, nBeta + idxMarge, stdNormal,
		      nBeta + idxMarge
		      FCONE FCONE FCONE FCONE);

      for (idxBeta=0;idxBeta<nBeta[idxMarge];idxBeta++)
	beta[cumBeta[idxMarge] + idxBeta] = stdNormal[idxBeta] +
	  conjMean[idxBeta];

      //The last step is to update the mean of the GP
      for (idxSite=0;idxSite<*nSite;idxSite++){
	GPmean[idxSite + idxMarge * *nSite] = 0;

	for (idxBeta=0;idxBeta<nBeta[idxMarge];idxBeta++)
	  GPmean[idxSite + idxMarge * *nSite] += dsgnMat[idxBeta * *nSite + idxSite +
							 cumBeta[idxMarge] * *nSite] *
	    beta[cumBeta[idxMarge] + idxBeta];
      }

      free(dummy);
      free(conjCovMat);
      free(conjCovMatChol);
      free(dummy2);
      free(conjMean);
      free(stdNormal);
    }

    /*----------------------------------------------------*/
    /*                                                    */
    /*        Updating the sills (conjugate prior)        */
    /*                                                    */
    /*----------------------------------------------------*/

    for (idxMarge=0;idxMarge<3;idxMarge++){
      if ((idxMarge==1) && *useLogLink){
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  resTop[idxSite] = log(gevParams[idxSite + idxMarge * *nSite]) -
	    GPmean[idxSite + idxMarge * *nSite];
      } else {
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  resTop[idxSite] = gevParams[idxSite + idxMarge * *nSite] -
	    GPmean[idxSite + idxMarge * *nSite];
      }

      // Compute resTop = covMatChol^(-T) %*% resTop
      F77_CALL(dtrsm)("L", "U", "T", "N", nSite, &oneInt, &one, covMatChol +
		      idxMarge * nSite2, nSite, resTop, nSite
		      FCONE FCONE FCONE FCONE);

      double shape = 0.5 * *nSite + hyperSill[2 * idxMarge];
      double scale = 0;
      for (idxSite=0;idxSite<*nSite;idxSite++)
	scale += resTop[idxSite] * resTop[idxSite];

      scale = hyperSill[1 + 2 * idxMarge] + 0.5 * sills[idxMarge] * scale;

      /* Rmk: If Y ~ Gamma(shape = shape, rate = 1 / scale) then X :=
	 1 / Y \sim IGamma(shape = shape, scale = scale) */
      sills[idxMarge] = 1 / rgamma(shape,  1 / scale);

      // Now we need to update the covariance matrix and its inverse
      switch(covmod[idxMarge]){
      case 1:
	flag = whittleMatern(distMat, nPairs, zero, sills[idxMarge], ranges[idxMarge],
			     smooths[idxMarge], covariances);
	break;
      case 2:
	flag = cauchy(distMat, nPairs, zero, sills[idxMarge], ranges[idxMarge],
		      smooths[idxMarge], covariances);
	break;
      case 3:
	flag = powerExp(distMat, nPairs, zero, sills[idxMarge], ranges[idxMarge],
			smooths[idxMarge], covariances);
	break;
      case 4:
	flag = bessel(distMat, nPairs, *dim, zero, sills[idxMarge], ranges[idxMarge],
		      smooths[idxMarge], covariances);
	break;
      }

      /* We need to fill in the upper triangular part of covMatChol with
	 covariances */
      {
	int current=-1;
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  for (idxSite2=idxSite;idxSite2<*nSite;idxSite2++){
	    current++;
	    covMatChol[idxSite + idxSite2 * *nSite + idxMarge * nSite2] = covariances[current];
	  }
      }

      // Cholesky decomposition of the covariance matrices
      F77_CALL(dpotrf)("U", nSite, covMatChol + idxMarge * nSite2, nSite,
		       &info FCONE);

      // Compute the log of the determinant of the proposal cov. mat.
      logDet[idxMarge] = 0;
      for (idxSite=0;idxSite<*nSite;idxSite++)
	logDet[idxMarge] += log(covMatChol[idxSite * (1 + *nSite) + idxMarge *
					   nSite2]);

      logDet[idxMarge] *= 2;
    }


    /*----------------------------------------------------*/
    /*                                                    */
    /*          Updating the ranges (M.-H. step)          */
    /*                                                    */
    /*----------------------------------------------------*/

    for (idxMarge=0;idxMarge<3;idxMarge++){
      if (propRanges[idxMarge] == 0)
	continue;

      double rangeProp = rlnorm(log(ranges[idxMarge]), propRanges[idxMarge]),
	logpropRatio = log(rangeProp / ranges[idxMarge]);

      switch(covmod[idxMarge]){
      case 1:
	flag = whittleMatern(distMat, nPairs, zero, sills[idxMarge], rangeProp,
			     smooths[idxMarge], covariances);
	break;
      case 2:
	flag = cauchy(distMat, nPairs, zero, sills[idxMarge], rangeProp,
		      smooths[idxMarge], covariances);
	break;
      case 3:
	flag = powerExp(distMat, nPairs, zero, sills[idxMarge], rangeProp,
			smooths[idxMarge], covariances);
	break;
      case 4:
	flag = bessel(distMat, nPairs, *dim, zero, sills[idxMarge], rangeProp,
		      smooths[idxMarge], covariances);
	break;
      }

      if (flag != 0){
	extRates[3 + idxMarge]++;
	continue;
      }

      /* We need to fill in the upper triangular part of covMatPropChol
	 with covariances */
      {
	int current=-1;
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  for (idxSite2=idxSite;idxSite2<*nSite;idxSite2++){
	    current++;
	    covMatPropChol[idxSite + idxSite2 * *nSite] = covariances[current];
	  }
      }

      // Cholesky decomposition of the proposal cov. mat.
      F77_CALL(dpotrf)("U", nSite, covMatPropChol, nSite, &info FCONE);

      if (info != 0){
	extRates[3 + idxMarge]++;
	continue;
      }

      // Log of the determinant of the proposal cov. mat.
      logDetProp = 0;
      for (idxSite=0;idxSite<*nSite;idxSite++)
	logDetProp += log(covMatPropChol[idxSite * (1 + *nSite)]);

      logDetProp *= 2;

      if ((idxMarge==1) && *useLogLink){
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  resBottom[idxSite] = log(gevParams[idxSite + idxMarge * *nSite]) -
	    GPmean[idxSite + idxMarge * *nSite];
      } else {
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  resBottom[idxSite] = gevParams[idxSite + idxMarge * *nSite] -
	    GPmean[idxSite + idxMarge * *nSite];
      }

      memcpy(resTop, resBottom, *nSite * sizeof(double));

      F77_CALL(dtrsm)("L", "U", "T", "N", nSite, &oneInt, &one, covMatChol +
		      idxMarge * nSite2, nSite, resBottom, nSite
		      FCONE FCONE FCONE FCONE);
      F77_CALL(dtrsm)("L", "U", "T", "N", nSite, &oneInt, &one, covMatPropChol,
		      nSite, resTop, nSite
		      FCONE FCONE FCONE FCONE);

      double top = logDetProp, bottom = logDet[idxMarge],
	logpriorRatio = (hyperRange[2 * idxMarge] - 1) *
	log(rangeProp / ranges[idxMarge]) + (ranges[idxMarge] - rangeProp) /
	hyperRange[2 * idxMarge + 1];

      for (idxSite=0;idxSite<*nSite;idxSite++){
	top += resTop[idxSite] * resTop[idxSite];
	bottom += resBottom[idxSite] * resBottom[idxSite];
      }

      top *= -0.5;
      bottom *= -0.5;

      if (unif_rand() < exp(top - bottom + logpriorRatio + logpropRatio)){
	ranges[idxMarge] = rangeProp;
	logDet[idxMarge] = logDetProp;
	memcpy(covMatChol + idxMarge * nSite2, covMatPropChol, nSite2 *
	       sizeof(double));
	accRates[3 + idxMarge]++;
      }
    }

    /*----------------------------------------------------*/
    /*                                                    */
    /*         Updating the smooths (M.-H. step)          */
    /*                                                    */
    /*----------------------------------------------------*/

    for (idxMarge=0;idxMarge<3;idxMarge++){
      if (propSmooths[idxMarge] == 0)
	continue;

      double smoothProp = rlnorm(log(smooths[idxMarge]), propSmooths[idxMarge]),
	logpropRatio = log(smoothProp / smooths[idxMarge]);

      switch(covmod[idxMarge]){
      case 1:
	flag = whittleMatern(distMat, nPairs, zero, sills[idxMarge], ranges[idxMarge],
			     smoothProp, covariances);
	break;
      case 2:
	flag = cauchy(distMat, nPairs, zero, sills[idxMarge], ranges[idxMarge],
		      smoothProp, covariances);
	break;
      case 3:
	flag = powerExp(distMat, nPairs, zero, sills[idxMarge], ranges[idxMarge],
			smoothProp, covariances);
	break;
      case 4:
	flag = bessel(distMat, nPairs, *dim, zero, sills[idxMarge], ranges[idxMarge],
		      smoothProp, covariances);
	break;
      }

      if (flag != 0){
    	extRates[6 + idxMarge]++;
    	continue;
      }

      /* We need to fill in the upper triangular part of covMatPropChol
    	 with covariances */
      {
    	int current=-1;
    	for (idxSite=0;idxSite<*nSite;idxSite++)
    	  for (idxSite2=idxSite;idxSite2<*nSite;idxSite2++){
    	    current++;
    	    covMatPropChol[idxSite + idxSite2 * *nSite] = covariances[current];
    	  }
      }

      // Cholesky decomposition of the proposal cov. mat.
      F77_CALL(dpotrf)("U", nSite, covMatPropChol, nSite, &info FCONE);

      if (info != 0){
    	extRates[6 + idxMarge]++;
    	continue;
      }

      // Log of the determinant of the proposal cov. mat.
      logDetProp = 0;
      for (idxSite=0;idxSite<*nSite;idxSite++)
    	logDetProp += log(covMatPropChol[idxSite * (1 + *nSite)]);

      logDetProp *= 2;

      if ((idxMarge==1) && *useLogLink){
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  resBottom[idxSite] = log(gevParams[idxSite + idxMarge * *nSite]) -
	    GPmean[idxSite + idxMarge * *nSite];
      } else {
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  resBottom[idxSite] = gevParams[idxSite + idxMarge * *nSite] -
	    GPmean[idxSite + idxMarge * *nSite];
      }

      memcpy(resTop, resBottom, *nSite * sizeof(double));

      F77_CALL(dtrsm)("L", "U", "T", "N", nSite, &oneInt, &one, covMatPropChol,
		      nSite, resTop, nSite
		      FCONE FCONE FCONE FCONE);
      F77_CALL(dtrsm)("L", "U", "T", "N", nSite, &oneInt, &one, covMatChol +
		      idxMarge * nSite2, nSite, resBottom, nSite
		      FCONE FCONE FCONE FCONE);

      double top = logDetProp, bottom = logDet[idxMarge],
    	logpriorRatio = (hyperSmooth[2 * idxMarge] - 1) *
    	log(smoothProp / smooths[idxMarge]) + (smooths[idxMarge] - smoothProp) /
    	hyperSmooth[2 * idxMarge + 1];

      for (idxSite=0;idxSite<*nSite;idxSite++){
    	top += resTop[idxSite] * resTop[idxSite];
    	bottom += resBottom[idxSite] * resBottom[idxSite];
      }

      top *= -0.5;
      bottom *= -0.5;

      if (unif_rand() < exp(top - bottom + logpriorRatio + logpropRatio)){
    	smooths[idxMarge] = smoothProp;
    	logDet[idxMarge] = logDetProp;
    	memcpy(covMatChol + idxMarge * nSite2, covMatPropChol, nSite2 *
    	       sizeof(double));
    	accRates[6 + idxMarge]++;
      }
    }

    iter++;

    //Need to store the new state into the mc object.
    if ((iter > *burnin) & ((iter % *thin) == 0)){
      mcLoc[nBeta[0] + iterThin * lagLoc] = sills[0];
      mcLoc[nBeta[0] + 1 + iterThin * lagLoc] = ranges[0];
      mcLoc[nBeta[0] + 2 + iterThin * lagLoc] = smooths[0];

      mcScale[nBeta[1] + iterThin * lagScale] = sills[1];
      mcScale[nBeta[1] + 1 + iterThin * lagScale] = ranges[1];
      mcScale[nBeta[1] + 2 + iterThin * lagScale] = smooths[1];

      mcShape[nBeta[2] + iterThin * lagShape] = sills[2];
      mcShape[nBeta[2] + 1 + iterThin * lagShape] = ranges[2];
      mcShape[nBeta[2] + 2 + iterThin * lagShape] = smooths[2];

      for (idxBeta=0;idxBeta<nBeta[0];idxBeta++)
	mcLoc[idxBeta + iterThin * lagLoc] = beta[idxBeta];

      for (idxBeta=0;idxBeta<nBeta[1];idxBeta++)
	mcScale[idxBeta + iterThin * lagScale] = beta[cumBeta[1] + idxBeta];

      for (idxBeta=0;idxBeta<nBeta[2];idxBeta++)
	mcShape[idxBeta + iterThin * lagShape] = beta[cumBeta[2] + idxBeta];

      for (idxSite=0;idxSite<*nSite;idxSite++){
	mcLoc[nBeta[0] + 3 + idxSite + iterThin * lagLoc] = gevParams[idxSite];
	mcScale[nBeta[1] + 3 + idxSite + iterThin * lagScale] = gevParams[*nSite + idxSite];
	mcShape[nBeta[2] + 3 + idxSite + iterThin * lagShape] = gevParams[2 * *nSite + idxSite];
      }
      iterThin++;
    }
  }
  PutRNGstate();

  for (int i=0;i<9;i++){
    accRates[i] /= (double) iter;
    extRates[i] /= (double) iter;
  }

  return;
}

void DIC(int *nChain, int *nSite, int *nObs, double *data, double *chainLoc,
	 double *chainScale, double *chainShape, double *postLoc,
	 double *postScale, double *postShape, double *dic, double *effNpar,
	 double *dbar){

  double tmp=0;
#pragma omp parallel for reduction(+:tmp)
  for (int i=0;i<*nChain;i++)
    for (int j=0;j<*nSite;j++){
      double dummy = 0;
      gevlik(data + j * *nObs, nObs, chainLoc + j * *nChain + i,
	     chainScale + j * *nChain + i, chainShape + j * *nChain + i,
	     &dummy);
      tmp += dummy;
    }

  *dbar = -2 * tmp / ((double) *nChain);

  tmp=0;
#pragma omp parallel for reduction(+:tmp)
  for (int i=0;i<*nSite;i++){
    double dummy = 0;
    gevlik(data + i * *nObs, nObs, postLoc + i, postScale + i, postShape + i,
	   &dummy);
    tmp += dummy;
  }

  *effNpar = *dbar + 2 * tmp;
  *dic = *effNpar + *dbar;

  return;
}

