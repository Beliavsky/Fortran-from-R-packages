#include "header.h"

double wlplikschlather(double *data, double *rho, double *jac,
		       int nObs, int nSite, double *weights){
  //This function computes the log-pairwise likelihood for the
  //schlather model.

  const int nPair = nSite * (nSite - 1) / 2;
  double dns = 0.0;

  for (int currentPair=0;currentPair<nPair;currentPair++){
    int i,j;
    getSiteIndex(currentPair, nSite, &i, &j);

    if (weights[currentPair] != 0){
      if (rho[currentPair] > .99999996)
	/* This means that only data1 or data2 contributes to the
	   log-likelihoood.  Hence we consider these parameters as
	   unfeasible as the bivariate distribution in this case
	   degenerates.

	   Rmq: a = .99999996 is the limiting numerical precision for
	   which a^2 = 1 */

	return (0.00000005 + rho[currentPair]) * (0.00000005 + rho[currentPair]) * MINF;


      for (int k=0;k<nObs;k++){
	if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	  continue;

	double data1Square = data[k + i * nObs] * data[k + i * nObs],
	  data2Square = data[k + j * nObs] * data[k + j * nObs],
	  twiceData1Data2 = 2 * data[k + i * nObs] * data[k + j * nObs],
	  c1 = sqrt(data1Square + data2Square - twiceData1Data2 * rho[currentPair]);

	//It's the log of the joint CDF
	double lFvec = - (data[k + i * nObs] + data[k + j * nObs] + c1) / twiceData1Data2;

	//It's the partial derivative for marge 1
	double  dvecM1 = -(rho[currentPair] * data[k + i * nObs] - c1 -
			   data[k + j * nObs]) / (2 * c1 * data1Square);

	//It's the partial derivative for marge 2
	double dvecM2 = -(rho[currentPair] * data[k + j * nObs] - c1 -
			  data[k + i * nObs]) / (2 * c1 * data2Square);

	//Rmq: to have dvecM1 and dvecM2 we have to multiply
	//them by Fvec[i]. It's not done yet as dvecMixed has to be
	//computed first.

	//It's the mixed partial derivative
	double dvecMixed = (1 - rho[currentPair] * rho[currentPair]) / (2 * c1 * c1 * c1) +
	  dvecM1 * dvecM2;

	//Now the final step, multiplying by Fvec and the gradient
	dns += weights[currentPair] * (log(dvecMixed) + lFvec + jac[k + i * nObs] +
				       jac[k + j * nObs]);
      }
    }
  }

  return dns;
}

double wlpliksmith(double *data, double *mahalDist, double *jac,
		   int nObs, int nSite, double *weights){
  //This function computes the log-pairwise likelihood for the
  //smith model.

  const int nPairs = nSite * (nSite - 1) / 2;
  double dns = 0.0;

  for (int currentPair=0;currentPair<nPairs;currentPair++){
    int i, j;
    getSiteIndex(currentPair, nSite, &i, &j);

    double imahal = 1 / mahalDist[currentPair];

    if (weights[currentPair] != 0){
      for (int k=0;k<nObs;k++){
	  if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	    continue;

	  double idata1 = 1 / data[k + i * nObs],
	    idata2 = 1 / data[k + j * nObs],
	    idata1Idata2Imahal = idata1 * idata2 * imahal,
	    c1 = log(data[k + j * nObs] * idata1) * imahal + 0.5 * mahalDist[currentPair],
	    c2 = mahalDist[currentPair] - c1;

	  if ((c1 > 38) && (c2 < -38))
	    // Contribution of site 1 only
	    dns += 2 * log(idata1) - idata1 + jac[k + i * nObs] + jac[k + j * nObs];

	  else if ((c1 < -38) && (c2 > 38))
	    // Contribution of site 2 only
	    dns += 2 * log(idata2) - idata2 + jac[k + i * nObs] + jac[k + j * nObs];

	  else if ((c1 > 38) && (c2 > 38))
	    // site 1 and site 2 are independent
	    dns += 2 * log(idata1 * idata2) - idata1 - idata2 + jac[k + i * nObs] + jac[k + j * nObs];

	  else {
	    // regular case (since the case c1 < -38 and c2 < -38 is impossible)

	    double dnormc1 = dnorm(c1, 0, 1, 0),
	      dnormc2 = dnorm(c2, 0, 1, 0),
	      pnormc1 = pnorm(c1, 0, 1, 1, 0),
	      pnormc2 = pnorm(c2, 0, 1, 1, 0),

	      //It's the log of the joint CDF
	      lFvec = -pnormc1 * idata1 - pnormc2 * idata2,

	      //It's the partial derivative for marge 1
	      dvecM1 = (dnormc1 * imahal + pnormc1) * idata1 * idata1 - dnormc2 *
	      idata1Idata2Imahal,

	      //It's the partial derivative for marge 2
	      dvecM2 = (dnormc2 * imahal + pnormc2) * idata2 * idata2 - dnormc1 *
	      idata1Idata2Imahal,

	      //Rmq: to have dvecM1 and dvecM2 we have to multiply
	      //them by Fvec[i]. It's not done yet as dvecMixed has to be
	      //computed first.

	      //It's the mixed partial derivative
	      dvecMixed = dvecM1 * dvecM2 + (data[k + j * nObs] * c2 * dnormc1 +
					     data[k + i * nObs] * c1 * dnormc2) *
	      idata1Idata2Imahal * idata1Idata2Imahal;

	  //Now the final step, multiplying by Fvec and the gradient
	    dns += weights[currentPair] * (log(dvecMixed) + lFvec + jac[k + i * nObs] +
					   jac[k + j * nObs]);
	}
      }
    }
  }

  return dns;
}

double wlplikschlatherind(double *data, double alpha, double *rho,
			  double *jac, int nObs, int nSite, double *weights){
  //This function computes the log-pairwise likelihood for the
  //schlather model allowing for independence.

  int i, j, k, currentPair = -1;
  double c1, dns = 0.0, lFvec, dvecM1, dvecM2, dvecMixed,
    data1Square, data2Square, twiceData1Data2;
  //c1 is a useful quantity - see documentation

  if (alpha == 0)
    //There is no independence part i.e. only a pure Schlather model
    dns = wlplikschlather(data, rho, jac, nObs, nSite, weights);

  else if (alpha == 1){
    //The process is a pure noise
    for (i=0;i<(nSite - 1);i++){
      for (j=i+1;j<nSite;j++){
	currentPair++;
	if (weights[currentPair] != 0){
	  for (k=nObs;k--;){
	    if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	      continue;
	    dns += weights[currentPair] * (- 1 / data[k + i * nObs] - 1 / data[k + j * nObs] -
					   2 * log(data[k + i * nObs] * data[k + j * nObs]) +
					   jac[k + i * nObs] + jac[k + j * nObs]);
	  }
	}
      }
    }
  }

  else {
    //This is a mixture between noise and Schlather
    for (i=0;i<(nSite-1);i++){
      for (j=i+1;j<nSite;j++){
	currentPair++;

	if (weights[currentPair] != 0){
	  if (rho[currentPair] > .99999996)
	    /* This means that only data1 or data2 contributes to the
	       log-likelihoood for the Schlather part.  Hence we
	       consider these parameters as unfeasible as the
	       bivariate distribution in this case degenerates.

	       Rmq: a = .99999996 is the limiting numerical precision for
	       which a^2 = 1 */

	    return (0.00000005 + rho[currentPair]) * (0.00000005 + rho[currentPair]) * MINF;


	  for (k=nObs;k--;){
	    if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	      continue;

	    data1Square = data[k + i * nObs] * data[k + i * nObs];
	    data2Square = data[k + j * nObs] * data[k + j * nObs];
	    twiceData1Data2 = 2 * data[k + i * nObs] * data[k + j * nObs];

	    c1 = sqrt(data1Square + data2Square - twiceData1Data2 * rho[currentPair]);

	    //It's the log of the joint CDF
	    lFvec = ((-alpha - 1) * (data[k + i * nObs] + data[k + j * nObs]) +
		     (alpha - 1) * c1) / twiceData1Data2;

	    //It's the partial derivative for marge 1
	    dvecM1 = (alpha - 1) * (rho[currentPair] * data[k + i * nObs] - c1 -
				    data[k + j * nObs]) / (2 * c1 * data1Square) +
	      alpha / data1Square;

	    //It's the partial derivative for marge 2
	    dvecM2 = (alpha - 1) * (rho[currentPair] * data[k + j * nObs] - c1 -
				    data[k + i * nObs]) / (2 * c1 * data2Square) +
	      alpha / data2Square;

	    //Rmq: to have dvecM1 and dvecM2 we have to multiply
	    //them by Fvec[i]. It's not done yet as dvecMixed has to be
	    //computed first.

	    //It's the mixed partial derivative
	    dvecMixed = (1 - alpha) * (1 - rho[currentPair] * rho[currentPair]) /
	      (2 * c1 * c1 * c1) + dvecM1 * dvecM2;

	    //Now the final step, multiplying by Fvec and the gradient
	    dns += weights[currentPair] * (log(dvecMixed) + lFvec + jac[k + i * nObs] +
					   jac[k + j * nObs]);

	  }
	}
      }
    }
  }

  return dns;
}

double wlplikextremalt(double *data, double *rho, double df, double *jac,
		       int nObs, int nSite, double *weights){
  //This function computes the log-pairwise likelihood for the
  //extremal t model.

  const int nPair = nSite * (nSite - 1) / 2;
  double dns = 0.0, idf = 1 / df, dfPlus1 = df + 1;

  for (int currentPair=0;currentPair<nPair;currentPair++){
    int i, j;
    getSiteIndex(currentPair, nSite, &i, &j);

      if (weights[currentPair] != 0){
	if (rho[currentPair] > .99999996)
	  /* Rmq: a = .99999996 is the limiting numerical precision for
	     which a^2 = 1 */
	  return (0.00000005 + rho[currentPair]) * (0.00000005 + rho[currentPair]) * MINF;

	double a = sqrt(dfPlus1 / (1 - rho[currentPair] * rho[currentPair]));

	for (int k=0;k<nObs;k++){
	  if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	    continue;

	  double idata1 = 1 / data[k + i * nObs],
	    idata2 = 1 / data[k + j * nObs],
	    data2_1 = R_pow(data[k + j * nObs] * idata1, idf),
	    data1_2 = 1 / data2_1,
	    c1 = (data2_1 - rho[currentPair]) * a,
	    c2 = (data1_2 - rho[currentPair]) * a,
	    dtc1 = dt(c1, dfPlus1, 0),
	    dtc2 = dt(c2, dfPlus1, 0);

	  if ((dtc1 == 0) || (dtc2 == 0))
	    //The bivariate distribution degenerates
	    return MINF;

	  double ptc1 = pt(c1, dfPlus1, 1, 0),
	    ptc2 = pt(c2, dfPlus1, 1, 0),
	  //It's the log of the joint CDF
	    lFvec = -ptc1 * idata1 - ptc2 * idata2,
	  //It's the partial derivative for marge 1
	  dvecM1 = idata1 * idata1 * ptc1 + idata1 * idata1 * a * idf *
	    data2_1 * dtc1 - idata1 * idata2 * a * idf * data1_2 * dtc2,
	  //It's the partial derivative for marge 2
	  dvecM2 = idata2 * idata2 * ptc2 + idata2 * idata2 * a * idf *
	    data1_2 * dtc2 - idata1 * idata2 * a * idf * data2_1 * dtc1;

	  //Rmq: to have dvecM1 and dvecM2 we have to multiply
	  //them by Fvec[i]. It's not done yet as dvecMixed has to be
	  //computed first.

	  //It's the mixed partial derivative
	  //Below are the derivative of the t density with df + 1 DoF at c1 and c2
	  double dertc1 = -(1 + 1 / dfPlus1) * c1 / (1  + c1 * c1 / dfPlus1) * dtc1,
	    dertc2 = -(1 + 1 / dfPlus1) * c2 / (1  + c2 * c2 / dfPlus1) * dtc2;
	  //Below is the mixed derivative of lFvec
	  double dlFvecMixed = idata1 * idata2 * idf * idf * a *
	    (dfPlus1 * (data2_1 * idata1 * dtc1 + data1_2 * idata2 * dtc2) +
	     (data2_1 * data2_1 * idata1 * dertc1 + data1_2 * data1_2 * idata2 * dertc2) *
	     a);

	  double dvecMixed = dlFvecMixed + dvecM1 * dvecM2;

	  //Now the final step, multiplying by Fvec and the gradient
	  dns += weights[currentPair] * (log(dvecMixed) + lFvec + jac[k + i * nObs] +
					 jac[k + j * nObs]);
	}
      }
    }

  return dns;
}
