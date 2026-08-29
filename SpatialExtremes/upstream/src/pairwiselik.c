#include "header.h"

double lplikschlather(double *data, double *rho, double *jac,
		      int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //schlather model.

  const int nPairs = nSite * (nSite - 1) / 2;
  double dns = 0.0;

#pragma omp parallel for reduction(+:dns)
  for (int currentPair=0;currentPair<nPairs;currentPair++){

    int i, j;
    getSiteIndex(currentPair, nSite, &i, &j);

    if (rho[currentPair] > .99999996){
      /* This means that only data1 or data2 contributes to the
	 log-likelihoood.

	 Rmq: a = .99999996 is the limiting numerical precision for
	 which a^2 = 1 */

      for (int k=0;k<nObs;k++){
	if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	  continue;

	if (data[k + i * nObs] >= data[k + j * nObs])
	  dns += -2 * log(data[k + j * nObs]) - 1 / data[k + j * nObs] + jac[k + i * nObs] + jac[k + j * nObs];

	else
	  dns += -2 * log(data[k + i * nObs]) - 1 / data[k + i * nObs] + jac[k + i * nObs] + jac[k + j * nObs];
      }
    }

    else {
      for (int k=0;k<nObs;k++){

	if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	  continue;

	double data1Square = data[k + i * nObs] * data[k + i * nObs],
	  data2Square = data[k + j * nObs] * data[k + j * nObs],
	  twiceData1Data2 = 2 * data[k + i * nObs] * data[k + j * nObs],
	  c1 = sqrt(data1Square + data2Square - twiceData1Data2 * rho[currentPair]),

	  //It's the log of the joint CDF
	  lFvec = - (data[k + i * nObs] + data[k + j * nObs] + c1) / twiceData1Data2,

	  //It's the partial derivative for marge 1
	  dvecM1 = -(rho[currentPair] * data[k + i * nObs] - c1 -
		     data[k + j * nObs]) / (2 * c1 * data1Square),

	  //It's the partial derivative for marge 2
	  dvecM2 = -(rho[currentPair] * data[k + j * nObs] - c1 -
		     data[k + i * nObs]) / (2 * c1 * data2Square),

	  //Rmq: to have dvecM1 and dvecM2 we have to multiply
	  //them by Fvec[i]. It's not done yet as dvecMixed has to be
	  //computed first.

	  //It's the mixed partial derivative
	  dvecMixed = (1 - rho[currentPair] * rho[currentPair]) / (2 * c1 * c1 * c1) +
	  dvecM1 * dvecM2;

	//Now the final step, multiplying by Fvec and the gradient
	dns += log(dvecMixed) + lFvec + jac[k + i * nObs] + jac[k + j * nObs];
      }
    }
  }

  return dns;
}

double lpliksmith(double *data, double *mahalDist, double *jac,
		  int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //smith model.

  const int nPairs = nSite * (nSite - 1) / 2;
  double dns = 0.0;

#pragma omp parallel for reduction(+:dns)
  for (int currentPair=0;currentPair<nPairs;currentPair++){
    int i, j;
    getSiteIndex(currentPair, nSite, &i, &j);

    double imahal = 1 / mahalDist[currentPair];

    for (int k=0;k<nObs;k++){

      if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	continue;

      double idata1 = 1 / data[k + i * nObs],
	idata2 = 1 / data[k + j * nObs],
	idata1Idata2Imahal = idata1 * idata2 * imahal,
	c1 = log(data[k + j * nObs] * idata1) * imahal + 0.5 * mahalDist[currentPair],
	c2 = mahalDist[currentPair] - c1;

      if ((c1 > 38) && (c2 < -38)){
	// Contribution of site 1 only
	dns += 2 * log(idata1) - idata1 + jac[k + i * nObs] + jac[k + j * nObs];
	//printf("case 1: mahal = %f\n", mahalDist[currentPair]);
      }

      else if ((c1 < -38) && (c2 > 38)){
	// Contribution of site 2 only
	dns += 2 * log(idata2) - idata2 + jac[k + i * nObs] + jac[k + j * nObs];
	//printf("case 2: mahal = %f\n", mahalDist[currentPair]);
      }

      else if ((c1 > 38) && (c2 > 38)){
	// site 1 and site 2 are independent
	dns += 2 * log(idata1 * idata2) - idata1 - idata2 + jac[k + i * nObs] + jac[k + j * nObs];
	//printf("case 3: mahal = %f\n", mahalDist[currentPair]);
      }

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
	dns += log(dvecMixed) + lFvec + jac[k + i * nObs] + jac[k + j * nObs];
      }
    }
  }

  return dns;
}

double lplikschlatherind(double *data, double alpha, double *rho,
			 double *jac, int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //schlather model allowing for independence.

  const int nPairs = nSite * (nSite - 1) / 2;
  double dns = 0.0;

  if (alpha == 0)
    //There is no independence part i.e. only a pure Schlather model
    dns = lplikschlather(data, rho, jac, nObs, nSite);

  else if (alpha == 1){
    //The process is a pure noise
#pragma omp parallel for reduction(+:dns)
    for (int currentPair=0;currentPair<nPairs;currentPair++){
      int i, j;
      getSiteIndex(currentPair, nSite, &i, &j);

      for (int k=0;k<nObs;k++){
	if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	  continue;

	dns += - 1 / data[k + i * nObs] - 1 / data[k + j * nObs] - 2 *
	  log(data[k + i * nObs] * data[k + j * nObs]) + jac[k + i * nObs] +
	  jac[k + j * nObs];
      }
    }
  }

  else {
    //This is a mixture between noise and Schlather
#pragma omp parallel for reduction(+:dns)
    for (int currentPair=0;currentPair<nPairs;currentPair++){
      int i, j;
      getSiteIndex(currentPair, nSite, &i, &j);

      if (rho[currentPair] > .99999996){
	/* This means that only data1 or data2 contributes to the
	   log-likelihoood for the Schlather part.

	   Rmq: a = .99999996 is the limiting numerical precision for
	   which a^2 = 1 */

	for (int k=0;k<nObs;k++){
	  if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	    continue;

	  if (data[k + i * nObs ] >= data[k + j * nObs])
	    dns += -2 * log(data[k + j * nObs]) - 1 / data[k + j * nObs] + jac[k + i * nObs] + jac[k + j * nObs];

	  else
	    dns += -2 * log(data[k + i * nObs]) - 1 / data[k + i * nObs] + jac[k + i * nObs] + jac[k + j * nObs];
	}
      }

      else {
	for (int k=0;k<nObs;k++){
	  if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	    continue;

	  double data1Square = data[k + i * nObs] * data[k + i * nObs],
	    data2Square = data[k + j * nObs] * data[k + j * nObs],
	    twiceData1Data2 = 2 * data[k + i * nObs] * data[k + j * nObs],
	    c1 = sqrt(data1Square + data2Square - twiceData1Data2 * rho[currentPair]),

	    //It's the log of the joint CDF
	    lFvec = ((-alpha - 1) * (data[k + i * nObs] + data[k + j * nObs]) +
		     (alpha - 1) * c1) / twiceData1Data2,

	    //It's the partial derivative for marge 1
	    dvecM1 = (alpha - 1) * (rho[currentPair] * data[k + i * nObs] - c1 -
				    data[k + j * nObs]) / (2 * c1 * data1Square) +
	    alpha / data1Square,

	    //It's the partial derivative for marge 2
	    dvecM2 = (alpha - 1) * (rho[currentPair] * data[k + j * nObs] - c1 -
				    data[k + i * nObs]) / (2 * c1 * data2Square) +
	    alpha / data2Square,

	    //Rmq: to have dvecM1 and dvecM2 we have to multiply
	    //them by Fvec[i]. It's not done yet as dvecMixed has to be
	    //computed first.

	    //It's the mixed partial derivative
	    dvecMixed = (1 - alpha) * (1 - rho[currentPair] * rho[currentPair]) /
	    (2 * c1 * c1 * c1) + dvecM1 * dvecM2;

	  //Now the final step, multiplying by Fvec and the gradient
	  dns += log(dvecMixed) + lFvec + jac[k + i * nObs] + jac[k + j * nObs];
	}
      }
    }
  }

  return dns;
}

double lplikextremalt(double *data, double *rho, double df, double *jac,
		      int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //extremal t model.

  const int nPairs = nSite * (nSite - 1) / 2;
  const double idf = 1 /df, dfPlus1 = df + 1;
  double dns = 0.0;

#pragma omp parallel for reduction(+:dns)
  for (int currentPair=0;currentPair<nPairs;currentPair++){
    int i, j;
    getSiteIndex(currentPair, nSite, &i, &j);

    if (rho[currentPair] > .99999996){
      /* Rmq: a = .99999996 is the limiting numerical precision for
	 which a^2 = 1 */

      for (int k=0;k<nObs;k++){
	if (ISNA(data[k + i * nObs]) || ISNA(data[k + j * nObs]))
	  continue;

	if (data[k + i * nObs ] >= data[k + j * nObs])
	  dns += -2 * log(data[k + j * nObs]) - 1 / data[k + j * nObs] + jac[k + i * nObs] + jac[k + j * nObs];

	else
	  dns += -2 * log(data[k + i * nObs]) - 1 / data[k + i * nObs] + jac[k + i * nObs] + jac[k + j * nObs];
      }
    }

    else {
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
	  dtc2 = dt(c2, dfPlus1, 0),
	  ptc1 = pt(c1, dfPlus1, 1, 0),
	  ptc2 = pt(c2, dfPlus1, 1, 0);

	if (ptc1 == 0)
	  //The bivariate distribution degenerates
	  dns +=  2 * log(idata2) - idata2 + jac[k + j * nObs];

	else if (ptc2 == 0)
	  //The bivariate distribution degenerates
	  dns += 2 * log(idata1) - idata1 + jac[k + i * nObs];

	else {
	  //It's the log of the joint CDF
	  double lFvec = -ptc1 * idata1 - ptc2 * idata2,

	    //It's the partial derivative for marge 1
	    dvecM1 = idata1 * (idata1 * ptc1 + a * idf *
			       (idata1 * data2_1 * dtc1 -
				idata2 * data1_2 * dtc2)),

	    //It's the partial derivative for marge 2
	    dvecM2 = idata2 * (idata2 * ptc2 + a * idf *
			       (idata2 * data1_2 * dtc2 -
				idata1 * data2_1 * dtc1)),

	    //Rmq: to have dvecM1 and dvecM2 we have to multiply
	    //them by Fvec[i]. It's not done yet as dvecMixed has to be
	    //computed first.

	    //It's the mixed partial derivative
	    //Below are the derivative of the t density with df + 1 DoF at c1 and c2
	    dertc1 = -(1 + 1 / dfPlus1) * c1 / (1  + c1 * c1 / dfPlus1) * dtc1,
	    dertc2 = -(1 + 1 / dfPlus1) * c2 / (1  + c2 * c2 / dfPlus1) * dtc2,

	    //Below is the mixed derivative of lFvec
	    dlFvecMixed = idata1 * idata2 * idf * idf * a *
	    (dfPlus1 * (data2_1 * idata1 * dtc1 + data1_2 * idata2 * dtc2) +
	     (data2_1 * data2_1 * idata1 * dertc1 + data1_2 * data1_2 * idata2 * dertc2) *
	     a),

	    dvecMixed = dlFvecMixed + dvecM1 * dvecM2;

	  //Now the final step, multiplying by Fvec and the gradient
	  dns += log(dvecMixed) + lFvec + jac[k + i * nObs] + jac[k + j * nObs];
	}
      }
    }
  }

  return dns;
}
