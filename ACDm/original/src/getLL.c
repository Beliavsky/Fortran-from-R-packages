#include "header.h" 
#include <Rinternals.h>
#include <Rdefines.h>

//START---getLL_dist-------------------------------------//
double getLL_dist(double *x,
		double *mu,
		double *e,
		int *N,
		int *dist,
		double *distPara,
		int *forceErrExpec){

	double LL = 0, gammaFncValue = 0, meanPara = 0;
	int i = 0;

	if(*forceErrExpec == 1){
		switch (*dist) {
		case 1: //Exponential distribution
			for (i = 0; i < *N; i++) {
				LL += -log(*(mu+i)) - *(e+i);
			}
			break;
		case 2: //Weibull distribution
			gammaFncValue = tgamma(1+1/(*distPara));
			for (i = 0; i < *N; i++) {
				LL += log(*distPara/(*(x+i))) + 
				*distPara*log(gammaFncValue*(*(x+i))/(*(mu+i))) - pow(gammaFncValue*(*(x+i))/(*(mu+i)), *distPara);
			}
			break;
		case 3: //Burr distribution
		
			meanPara = lgamma(1 + 1 / distPara[0]);
			meanPara += lgamma(1 / distPara[1] - 1 / distPara[0]);
			meanPara -= lgamma(1 + 1 / distPara[1]);
			meanPara = exp(meanPara);
			meanPara /= pow(distPara[1], (1 + 1 / distPara[0]));
			meanPara = pow(meanPara, distPara[0]);

			LL = (log(meanPara * distPara[0])) * (*N); //multiplies the constant terms with N outside of the loop below
			for (i = 0; i < *N; i++) {
				LL += -log(mu[i] * e[i]) + distPara[0] * log(e[i]) - 
				(1 + 1 / distPara[1]) * log(1 + distPara[1] * meanPara * pow(e[i], distPara[0]));
			}
			break;
		case 4: //Generalized Gamma distribution
		{
			double theta,  lnDelta, alphDelt, lnThet, lnThetGamm, thetaPow;

			theta = exp(lgamma(*distPara) - lgamma(*distPara + 1/(*(distPara + 1)))); //scale parameter
			lnDelta = log(*(distPara + 1));
			alphDelt = (*distPara) * (*(distPara + 1)) - 1;
			lnThet = -alphDelt * log(theta);
			lnThetGamm = log(theta) + lgamma(*distPara);
			thetaPow = pow(theta, -(*(distPara + 1)));
			LL = (lnDelta + lnThet - lnThetGamm) * (*N); //multiplies the constant terms with N outside of the loop below
			for (i = 0; i < *N; i++) {
				LL += alphDelt * log(*(e+i)) - pow((*(e+i)), *(distPara + 1)) * thetaPow - log(*(mu+i));
			}
			break;
		}
		case 5: //Generalized F distribution
		{
			double lnbetaFnPart, am, amLess1, etaPm, zetaInv;
			lnbetaFnPart = lgamma(*distPara + (*(distPara+1))) - lgamma(*distPara) - lgamma(*(distPara+1));
			am = *(distPara + 2) * (*distPara);
			amLess1 = am - 1;
			etaPm = *(distPara + 1) + (*distPara);
			zetaInv = pow(*(distPara + 1), -1/(*(distPara + 2)))
															* exp(lgamma(*distPara) + lgamma(*(distPara + 1)) - lgamma(*distPara + 1/(*(distPara + 2))) - lgamma(*(distPara + 1) - 1/(*(distPara + 2))));
			LL = (lnbetaFnPart + log(*(distPara + 2)) - am*log(zetaInv) + *(distPara + 1)*log(*(distPara + 1))) * (*N); //multiplies the constant terms with N outside of the loop below
			for (i = 0; i < *N; i++) {
				LL += -am * log(*(mu+i)) + amLess1 * log(*(x+i)) - etaPm * log(*(distPara + 1) + pow(*(e+i)/zetaInv, *(distPara + 2)));
			}
			break;
		}
		case 6: //q-Weibull
		{
			double a = *distPara, q = *(distPara + 1), b = 0, aLess1 = 0;
			b = lgamma(1/(q-1)) - lgamma(1/a) - lgamma(1/(q-1)-1/a-1);
			b = exp(b) * a * pow(q-1, (1+a)/a) / (2 - q);
			aLess1 = a - 1;

			LL = log(a*(2-q)/pow(b, a)) * (*N);
			for (i = 0; i < *N; i++) {
				LL += aLess1 * log(*(e+i)) + (1/(1-q)) * log(1 - (1-q) * pow(*(e+i)/b, a)) - log(*(mu+i));
			}
			break;
		}
		case 7: //q-Weibull and exponential discreetly mixed
		{
			double p = *distPara, a = *(distPara + 1), q = *(distPara + 2), lambda = *(distPara + 3), b = 0, const1 = 0, const2 = 0;
			b = lgamma(1/(q-1)) - lgamma(1/a) - lgamma(1/(q-1)-1/a-1);
			b = exp(b) * a * pow(q-1, (1+a)/a) / (2 - q);
			b = b * (1 - (1 - p) * lambda) / p;

			const1 = p * (2 - q) * a * pow(b, -a);
			const2 = (1 - p) / lambda;

			for (i = 0; i < *N; i++) {
				LL += log( const1 * pow(e[i], a-1) * pow(1 - (1 - q) * pow(e[i]/b, a), 1/(1-q)) + const2 * exp(-e[i]/lambda) ) - log(mu[i]);
			}

			break;
		}
		case 8: //q-Weibull and Weibull discreetly mixed
		{
			double p = *distPara, a = *(distPara + 1), q = *(distPara + 2), theta = *(distPara + 3), gamma = *(distPara + 4), b = 0, const1 = 0, const2 = 0;
			b = lgamma(1/(q-1)) - lgamma(1/a) - lgamma(1/(q-1)-1/a-1);
			b = exp(b) * a * pow(q-1, (1+a)/a) / (2 - q);
			b = b * (1 - (1 - p) * pow(theta, -1/gamma) * tgamma(1/gamma + 1)) / p;

			const1 = p * (2 - q) * a * pow(b, -a);
			const2 = (1 - p) * theta * gamma;

			for (i = 0; i < *N; i++) {
				LL += log( const1 * pow(e[i], a-1) * pow(1 - (1 - q) * pow(e[i]/b, a), 1/(1-q)) + const2 * pow(e[i], gamma - 1) * exp(-theta * pow(e[i], gamma) ) ) - log(mu[i]);
			}

			break;
		}
		case 9: //"mixinvgauss" - finite inverse Gaussian mixature
		{
			double theta = *distPara, lambda = *(distPara + 1), gamma = *(distPara + 2); //, const1 = 0, const2 = 0, const3 = 0;
			double phi = theta * (1 + theta * theta / lambda / (1 + gamma));
			//const1 = lambda/2;
			//const2 = gamma/phi;
			//const3 = phi * theta * theta;

			//LL = (log(phi*phi/(gamma + theta)) - .5 * (log(2 * M_PI / lambda)) - 1.5 * log(phi)) * (*N); //M_PI is pi
			LL=0;
			for (i = 0; i < *N; i++) {
				//LL += log(const2 + e[i]) - 1.5 * log(e[i]) - const1 * (phi * e[i] - theta) * (phi * e[i] - theta) / (e[i] * const3) - log(mu[i]);
				LL += -log(mu[i]) + log(gamma * mu[i] + phi * x[i]) - log(gamma + theta) + .5 * (log(lambda) + log(mu[i]) - log(2 * M_PI) - log(phi)) - 1.5 * log(x[i]) - 0.5 * lambda * (phi*x[i]-theta*mu[i]) * (phi*x[i]-theta*mu[i]) /(phi * x[i]*theta*theta*mu[i]);
			}
			break;
		}
		case 10: //Birnbaum Saunders
		{
			double kappa = *distPara;
			double const1 = 0;
					
			const1 = 1/(2 * kappa * kappa);
			
			LL = -log(2 * sqrt(2 * M_PI) * kappa) * (*N); //M_PI is pi
			for (i = 0; i < *N; i++) {
				//LL += -log(mu[i]) + log(pow(e[i], -1/2) + pow(e[i], -3/2)) - const1 * (e[i] + 1/e[i] - 2);
				LL += -log(mu[i]) + log(pow(mu[i]/x[i], 1/2) + pow(mu[i]/x[i], 3/2)) - const1 * (x[i]/mu[i] + mu[i]/x[i] - 2);
			}
			break;
		}
		default:
			break;
		}
	} else if(*forceErrExpec == 0){ //mean parameter is fixed - mean is not forced to be 1
		switch (*dist) {
		case 1: //Exponential distribution
			for (i = 0; i < *N; i++) {
				LL += -log(*(mu+i)) - *(e+i);
			}
			break;
		case 2: //Weibull distribution
			LL = log(distPara[0])  * (*N); //multiplies the constant terms with N outside of the loop below
			for (i = 0; i < *N; i++) {
				LL += -log(mu[i]) + (distPara[0] - 1) * log(e[i]) - pow(e[i], distPara[0]);
			}
			break;
		case 3: //Burr distribution
			meanPara = 1;

			LL = (log(meanPara * distPara[0])) * (*N); //multiplies the constant terms with N outside of the loop below
			for (i = 0; i < *N; i++) {
				LL += -log(mu[i]) + (distPara[0] - 1) * log(e[i]) - (1 + 1 / distPara[1]) * log(1 + distPara[1] * meanPara * pow(e[i], distPara[0]));
			}
			break;
		case 4: //Generalized Gamma distribution
		{
			double lnDelta, alphDelt, lnGamm;
			lnDelta = log(*(distPara + 1));
			alphDelt = (*distPara) * (*(distPara + 1)) - 1;
			lnGamm = lgamma(*distPara);
			LL = (lnDelta - lnGamm) * (*N); //multiplies the constant terms with N outside of the loop below
			for (i = 0; i < *N; i++) {
				LL += alphDelt * log(*(e+i)) - pow((*(e+i)), *(distPara + 1)) - log(*(mu+i));
			}
			break;
		}
		case 5: //Generalized F distribution
		{
			double lnbetaFnPart, am, amLess1, etaPm;
			lnbetaFnPart = lgamma(*distPara + (*(distPara+1))) - lgamma(*distPara) - lgamma(*(distPara+1));
			am = *(distPara + 2) * (*distPara);
			amLess1 = am - 1;
			etaPm = *(distPara + 1) + (*distPara);
			LL = (lnbetaFnPart + log(*(distPara + 2)) + *(distPara + 1)*log(*(distPara + 1))) * (*N); //multiplies the constant terms with N outside of the loop below
			for (i = 0; i < *N; i++) {
				LL += -am * log(*(mu+i)) + amLess1 * log(*(x+i)) - etaPm * log(*(distPara + 1) + pow(*(e+i), *(distPara + 2)));
			}
			break;
		}
		case 6: //q-Weibull
				{
					double a = *distPara, q = *(distPara + 1), aLess1 = 0;
					aLess1 = a - 1;

					LL = log(a * (2 - q)) * (*N);
					for (i = 0; i < *N; i++) {
						LL += aLess1 * log(*(e+i)) + (1/(1-q)) * log(1 - (1-q) * pow(*(e+i), a)) - log(*(mu+i));
					}
					break;
				}
		case 7: //q-Weibull and exponential discreetly mixed
		{
			double p = *distPara, a = *(distPara + 1), q = *(distPara + 2), lambda = *(distPara + 3), b = 0, const1 = 0, const2 = 0;
			b = 1;

			const1 = p*(2-q)*a/pow(b, a);
			const2 = (1-p)/lambda;

			for (i = 0; i < *N; i++) {
				LL += log( const1 * pow(e[i], a-1) * pow(1 - (1 - q) * pow(e[i]/b, a), 1/(1-q)) + const2 * exp(-e[i]/lambda) ) - log(mu[i]);
			}

			break;
		}
		case 8: //q-Weibull and Weibull discreetly mixed
		{
			double p = *distPara, a = *(distPara + 1), q = *(distPara + 2), theta = *(distPara + 3), gamma = *(distPara + 4), b = 0, const1 = 0, const2 = 0;
			b = 1;
			const2 = (1 - p) * theta * gamma;

			for (i = 0; i < *N; i++) {
				LL += log( const1 * pow(e[i], a-1) * pow(1 - (1 - q) * pow(e[i]/b, a), 1/(1-q)) + const2 * pow(e[i], gamma - 1) * exp(-theta * pow(e[i], gamma) ) ) - log(mu[i]);
			}

			break;
		}
		case 9: //"mixinvgauss" - finite inverse Gaussian mixature
		{
			double theta = *distPara, lambda = *(distPara + 1), gamma = *(distPara + 2); //, const1 = 0, const2 = 0, const3 = 0;
			double phi = 1;
			//const1 = lambda/2;
			//const2 = gamma/phi;
			//const3 = phi * theta * theta;

			//LL = (log(phi*phi/(gamma + theta)) - .5 * (log(2 * M_PI / lambda)) - 1.5 * log(phi)) * (*N); //M_PI is pi
			LL=0;
			for (i = 0; i < *N; i++) {
				//LL += log(const2 + e[i]) - 1.5 * log(e[i]) - const1 * (phi * e[i] - theta) * (phi * e[i] - theta) / (e[i] * const3) - log(mu[i]);
				LL += -log(mu[i]) + log(gamma * mu[i] + phi * x[i]) - log(gamma + theta) + .5 * (log(lambda) + log(mu[i]) - log(2 * M_PI) - log(phi)) - 1.5 * log(x[i]) - 0.5 * lambda * (phi*x[i]-theta*mu[i]) * (phi*x[i]-theta*mu[i]) /(phi * x[i]*theta*theta*mu[i]);
			}
			break;

		}
		case 10: //Birnbaum Saunders
		{
			double kappa = *distPara;
			double const1 = 0;

			const1 = 1/(2 * kappa * kappa);

			LL = -log(2 * sqrt(2 * M_PI) * kappa) * (*N); //M_PI is pi
			for (i = 0; i < *N; i++) {
				LL += -log(mu[i]) + log(pow(e[i], -1/2) + pow(e[i], -3/2)) - const1 * (e[i] + 1/e[i] - 2);
			}
			break;
		}
		default:
			break;
		}
	}
	return(LL);
}
//END---getLL_dist-------------------------------------//

//START---getLL_ACDcall----------------------------//
SEXP getLL_ACDcall(
		SEXP x,
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){


	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	int maxpq = max(p, q);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpq;

	double *px;
	int *pnewDay;
	px = REAL(x);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(mean)[0];
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(par)[0]; //adds the constant
			for(j = 1; j <= p; j++) pmu[i] += REAL(par)[j] * px[i - j]; //adds the lagged durations
			for(j = 1; j <= q; j++) pmu[i] += REAL(par)[j + p] * pmu[i - j]; //adds the lagged mus
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, N);
	} while (stopIndex != N);
	
	
	
	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_ACDcall----------------------------//

//START---getLL_ACDcallEx----------------------------//
SEXP getLL_ACDcallEx(
		SEXP x,
		SEXP z, //the exogenous regressors
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){


	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	int maxpq = max(p, q);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpq;

	double *px, *pz;
	int *pnewDay;
	px = REAL(x);
	pz = REAL(z);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	int k = (int)(length(z)/N); //the number of exogenous variables 
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(mean)[0];
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(par)[0]; //adds the constant
			for(j = 1; j <= p; j++) pmu[i] += REAL(par)[j] * px[i - j]; //adds the lagged durations
			for(j = 1; j <= q; j++) pmu[i] += REAL(par)[j + p] * pmu[i - j]; //adds the lagged mus
			for(j = 0; j <= k - 1; j++) pmu[i] += REAL(par)[j + 1 + q + p] * pz[i + j * N]; //adds the exogenous variables
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, N);
	} while (stopIndex != N);
	
	
	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_ACDcallEx----------------------------//

//START---getLL_EXACDcall----------------------------//
SEXP getLL_EXACDcall(
    SEXP x,
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], q = INTEGER(order)[1];
  int maxpq = max(p, q);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpq;
  
  double *px;
  int *pnewDay;
  px = REAL(x);
  pnewDay = INTEGER(newDay);
  
  int N = length(x), NnewDays = length(newDay);
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  double logMu[N];
  
  do {
    //at the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
    for (i = startIndex; i < stopIndex; i++){
      logMu[i] = log(REAL(mean)[0]);
      pmu[i] = exp(logMu[i]);
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    for(i = startIndex; i < stopIndex; i++){
      logMu[i] = REAL(par)[0]; //adds the constant
      for(j = 1; j <= p; j++) logMu[i] += REAL(par)[j] * presi[i - j] + REAL(par)[j + p] * fabs(presi[i - j] - 1); //adds the lagged durations
      for(j = 1; j <= q; j++) logMu[i] += REAL(par)[j + 2 *p] * logMu[i - j]; //adds the lagged mus
      pmu[i] = exp(logMu[i]);
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpq, N);
  } while (stopIndex != N);
  
  
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_EXACDcall----------------------------//

//START---getLL_EXACDcallEx----------------------------//
SEXP getLL_EXACDcallEx(
    SEXP x,
    SEXP z, //the exogenous regressors
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], q = INTEGER(order)[1];
  int maxpq = max(p, q);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpq;
  
  double *px, *pz;
  int *pnewDay;
  px = REAL(x);
  pz = REAL(z);
  pnewDay = INTEGER(newDay);
  
  int N = length(x), NnewDays = length(newDay);
  int k = (int)(length(z)/N); //the number of exogenous variables 
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  double logMu[N];
  
  do {
    //at the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
    for (i = startIndex; i < stopIndex; i++){
      logMu[i] = log(REAL(mean)[0]);
      pmu[i] = exp(logMu[i]);
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    for(i = startIndex; i < stopIndex; i++){
      logMu[i] = REAL(par)[0]; //adds the constant
      for(j = 1; j <= p; j++) logMu[i] += REAL(par)[j] * presi[i - j] + REAL(par)[j + p] * fabs(presi[i - j] - 1); //adds the lagged durations
      for(j = 1; j <= q; j++) logMu[i] += REAL(par)[j + 2 * p] * logMu[i - j]; //adds the lagged mus
      for(j = 0; j <= k - 1; j++) logMu[i] += REAL(par)[j + 1 + q + 2 * p] * pz[i + j * N]; //adds the exogenous variables
      pmu[i] = exp(logMu[i]);
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpq, N);
  } while (stopIndex != N);
  
  
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_EXACDcallEx----------------------------//

//START---getLL_LACD1call----------------------------//
SEXP getLL_LACD1call(
    SEXP x,
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], q = INTEGER(order)[1];
  int maxpq = max(p, q);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpq;
  
  double *px;
  int *pnewDay;
  px = REAL(x);
  pnewDay = INTEGER(newDay);
  
  int N = length(x), NnewDays = length(newDay);
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  double logMu[N];
  
  do {
    //in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
    for (i = startIndex; i < stopIndex; i++){
      logMu[i] = log(REAL(mean)[0]);
      pmu[i] = exp(logMu[i]);
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    for(i = startIndex; i < stopIndex; i++){
      logMu[i] = REAL(par)[0]; //adds the constant
      for(j = 1; j <= p; j++) logMu[i] += REAL(par)[j] * log(presi[i - j]); //adds the lagged durations
      for(j = 1; j <= q; j++) logMu[i] += REAL(par)[j + p] * logMu[i - j]; //adds the lagged mus
      pmu[i] = exp(logMu[i]);
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpq, N);
  } while (stopIndex != N);
  
  
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_LACD1call----------------------------//

//START---getLL_LACD1callEx----------------------------//
SEXP getLL_LACD1callEx(
		SEXP x,
		SEXP z, //the exogenous regressors
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){


	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	int maxpq = max(p, q);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpq;

	double *px, *pz;
	int *pnewDay;
	px = REAL(x);
	pz = REAL(z);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	int k = (int)(length(z)/N); //the number of exogenous variables 
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
	double logMu[N];

	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			logMu[i] = log(REAL(mean)[0]);
			pmu[i] = exp(logMu[i]);
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			logMu[i] = REAL(par)[0]; //adds the constant
			for(j = 1; j <= p; j++) logMu[i] += REAL(par)[j] * log(presi[i - j]); //adds the lagged durations
			for(j = 1; j <= q; j++) logMu[i] += REAL(par)[j + p] * logMu[i - j]; //adds the lagged mus
			for(j = 0; j <= k - 1; j++) logMu[i] += REAL(par)[j + 1 + q + p] * pz[i + j * N]; //adds the exogenous variables			
			pmu[i] = exp(logMu[i]);
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, N);
	} while (stopIndex != N);
	
	
	
	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_LACD1callEx----------------------------//

//START---getLL_BCACDcall----------------------------//
SEXP getLL_BCACDcall(
    SEXP x,
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], q = INTEGER(order)[1];
  int maxpq = max(p, q);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpq;
  
  double *px;
  int *pnewDay;
  px = REAL(x);
  pnewDay = INTEGER(newDay);
  
  int N = length(x), NnewDays = length(newDay);
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  
  double LogMu[N]; 
  double d = REAL(par)[1 + p + q];
  double a[p], b[q];
  
  for(j = 0; j < p; j++)
    a[j] = REAL(par)[j + 1];	
  for(j = 0; j < q; j++)
    b[j] = REAL(par)[j + 1 + p];	
  
  do {
    //in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
    for (i = startIndex; i < stopIndex; i++){
      pmu[i] = REAL(mean)[0];
      presi[i] = px[i]/pmu[i];
      LogMu[i] = log(REAL(mean)[0]);
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    for(i = startIndex; i < stopIndex; i++){
      LogMu[i] = REAL(par)[0]; //adds the constant
      for(j = 0; j < p; j++) LogMu[i] += a[j] * pow(presi[i - j - 1], d); //adds the p-part
      for(j = 0; j < q; j++) LogMu[i] += b[j] * LogMu[i - j - 1]; //adds the q-part				
      
      pmu[i] = exp(LogMu[i]);
      presi[i] = px[i]/pmu[i];
    }
    
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpq, N);
  } while (stopIndex != N);
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_BCACDcall----------------------------//


//START---getLL_BCACDcallEx----------------------------//
SEXP getLL_BCACDcallEx(
    SEXP x,
    SEXP z, //the exogenous regressors
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], q = INTEGER(order)[1];
  int maxpq = max(p, q);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpq;
  
  double *px, *pz;
  int *pnewDay;
  px = REAL(x);
  pz = REAL(z);
  pnewDay = INTEGER(newDay);
  
  int N = length(x), NnewDays = length(newDay);
  int k = (int)(length(z)/N); //the number of exogenous variables 
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  
  double LogMu[N]; //the mus to the d[0] power
  double d = REAL(par)[1 + p + q];
  double a[p], b[q];
  
  for(j = 0; j < p; j++)
    a[j] = REAL(par)[j + 1];	
  for(j = 0; j < q; j++)
    b[j] = REAL(par)[j + 1 + p];	
  
  do {
    //in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
    for (i = startIndex; i < stopIndex; i++){
      pmu[i] = REAL(mean)[0];
      presi[i] = px[i]/pmu[i];
      LogMu[i] = log(REAL(mean)[0]);
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    for(i = startIndex; i < stopIndex; i++){
      LogMu[i] = REAL(par)[0]; //adds the constant
      for(j = 0; j < p; j++) LogMu[i] += a[j] * pow(presi[i - j - 1], d); //adds the p-part
      for(j = 0; j < q; j++) LogMu[i] += b[j] * LogMu[i - j - 1]; //adds the q-part				
      for(j = 0; j <= k - 1; j++) LogMu[i] += REAL(par)[j + 1 + q + p] * pz[i + j * N]; //adds the exogenous variables		
      
      pmu[i] = exp(LogMu[i]);
      presi[i] = px[i]/pmu[i];
    }
    
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpq, N);
  } while (stopIndex != N);
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_BCACDcallEx----------------------------//

//START---getLL_BACDcall----------------------------//
SEXP getLL_BACDcall(
		SEXP x,
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){


	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	int maxpq = max(p, q);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpq;

	double *px;
	int *pnewDay;
	px = REAL(x);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

	double poweredMu[N]; //the mus to the d[0] power
	double d[2] = {REAL(par)[1 + p + q], REAL(par)[2 + p + q]};
	double a[p], b[q];

	for(j = 0; j < p; j++)
		a[j] = REAL(par)[j + 1];	
	for(j = 0; j < q; j++)
		b[j] = REAL(par)[j + 1 + p];	
	
	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(mean)[0];
			presi[i] = px[i]/pmu[i];
			poweredMu[i] = pow(REAL(mean)[0], d[0]);
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			poweredMu[i] = REAL(par)[0]; //adds the constant
			for(j = 0; j < p; j++) poweredMu[i] += a[j] * pow(presi[i - j - 1], d[1]); //adds the p-part
			for(j = 0; j < q; j++) poweredMu[i] += b[j] * poweredMu[i - j - 1]; //adds the q-part				

			pmu[i] = pow(poweredMu[i], 1/d[0]);
			presi[i] = px[i]/pmu[i];
		}
					
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, N);
	} while (stopIndex != N);

	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_BACDcall----------------------------//

//START---getLL_BACDcallEx----------------------------//
SEXP getLL_BACDcallEx(
		SEXP x,
		SEXP z, //the exogenous regressors
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){


	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	int maxpq = max(p, q);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpq;

	double *px, *pz;
	int *pnewDay;
	px = REAL(x);
	pz = REAL(z);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	int k = (int)(length(z)/N); //the number of exogenous variables 
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

	double poweredMu[N]; //the mus to the d[0] power
	double d[2] = {REAL(par)[1 + p + q], REAL(par)[2 + p + q]};
	double a[p], b[q];

	for(j = 0; j < p; j++)
		a[j] = REAL(par)[j + 1];	
	for(j = 0; j < q; j++)
		b[j] = REAL(par)[j + 1 + p];	
	
	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(mean)[0];
			presi[i] = px[i]/pmu[i];
			poweredMu[i] = pow(REAL(mean)[0], d[0]);
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			poweredMu[i] = REAL(par)[0]; //adds the constant
			for(j = 0; j < p; j++) poweredMu[i] += a[j] * pow(presi[i - j - 1], d[1]); //adds the p-part
			for(j = 0; j < q; j++) poweredMu[i] += b[j] * poweredMu[i - j - 1]; //adds the q-part				
			for(j = 0; j <= k - 1; j++) poweredMu[i] += REAL(par)[j + 1 + q + p] * pz[i + j * N]; //adds the exogenous variables			
			
			pmu[i] = pow(poweredMu[i], 1/d[0]);
			presi[i] = px[i]/pmu[i];
		}
					
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, N);
	} while (stopIndex != N);

	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_BACDcallEx----------------------------//

//START---getLL_ABACDcall----------------------------//
SEXP getLL_ABACDcall(
		SEXP x,
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){

	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	int maxpq = max(p, q);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpq;

	double *px;
	int *pnewDay;
	px = REAL(x);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

	double poweredMu[N]; //the mus to the d[0] power
	
	double c = REAL(par)[1+p+q];
	double v = REAL(par)[2+p+q];
	double d[2] = {REAL(par)[3+p+q], REAL(par)[4+p+q]};	
	
	double a[p], b[q];

	for(j = 0; j < p; j++){
		a[j] = REAL(par)[j + 1];
	}
	for(j = 0; j < q; j++)
		b[j] = REAL(par)[j + 1 +  p];

	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(mean)[0];
			presi[i] = 1;
			poweredMu[i] = pow(REAL(mean)[0], d[0]);
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			poweredMu[i] = REAL(par)[0]; //adds the constant
			for(j = 0; j < p; j++) poweredMu[i] += a[j]*pow(fabs(presi[i - j - 1] - v) + c * (presi[i - j - 1] - v), d[1]); //adds the p-part
			for(j = 0; j < q; j++) poweredMu[i] += b[j]*poweredMu[i - j - 1]; //adds the q-part
			
			pmu[i] = pow(poweredMu[i], 1/d[0]);
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, N);
	} while (stopIndex != N);

	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_ABACDcall----------------------------//

//START---getLL_ABACDcallEx----------------------------//
SEXP getLL_ABACDcallEx(
		SEXP x,
		SEXP z, //the exogenous regressors
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){

	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	int maxpq = max(p, q);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpq;

	double *px, *pz;
	int *pnewDay;
	px = REAL(x);
	pz = REAL(z);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	int k = (int)(length(z)/N); //the number of exogenous variables 
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

	double poweredMu[N]; //the mus to the d[0] power
	
	double c = REAL(par)[1+p+q];
	double v = REAL(par)[2+p+q];
	double d[2] = {REAL(par)[3+p+q], REAL(par)[4+p+q]};	
	
	double a[p], b[q];

	for(j = 0; j < p; j++){
		a[j] = REAL(par)[j + 1];
	}
	for(j = 0; j < q; j++)
		b[j] = REAL(par)[j + 1 +  p];

	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(mean)[0];
			presi[i] = 1;
			poweredMu[i] = pow(REAL(mean)[0], d[0]);
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			poweredMu[i] = REAL(par)[0]; //adds the constant
			for(j = 0; j < p; j++) poweredMu[i] += a[j]*pow(fabs(presi[i - j - 1] - v) + c * (presi[i - j - 1] - v), d[1]); //adds the p-part
			for(j = 0; j < q; j++) poweredMu[i] += b[j]*poweredMu[i - j - 1]; //adds the q-part
			for(j = 0; j <= k - 1; j++) poweredMu[i] += REAL(par)[j + 1 + q + p] * pz[i + j * N]; //adds the exogenous variables			
			
			pmu[i] = pow(poweredMu[i], 1/d[0]);
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, N);
	} while (stopIndex != N);

	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_ABACDcallEx----------------------------//

//START---getLL_AACDcall----------------------------//
SEXP getLL_AACDcall(
    SEXP x,
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], q = INTEGER(order)[1];
  int maxpq = max(p, q);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpq;
  
  double *px;
  int *pnewDay;
  px = REAL(x);
  pnewDay = INTEGER(newDay);
  
  int N = length(x), NnewDays = length(newDay);
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  
  double poweredMu[N]; //the mus to the d[0] power
  double c = REAL(par)[1 + p + q];
  double v = REAL(par)[2 + p + q];
  double d[2] = {REAL(par)[p + q + 3], REAL(par)[p + q + 4]};
  double a[p], b[q];
  
  for(j = 0; j < p; j++) a[j] = REAL(par)[j + 1];
  for(j = 0; j < q; j++) b[j] = REAL(par)[j + 1 + p];
  
  do {
    //in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
    for (i = startIndex; i < stopIndex; i++){
      pmu[i] = REAL(mean)[0];
      presi[i] = 1;
      poweredMu[i] = pow(REAL(mean)[0], d[0]);
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    for(i = startIndex; i < stopIndex; i++){
      poweredMu[i] = REAL(par)[0]; //adds the constant
      for(j = 0; j < p; j++) poweredMu[i] += a[j] * poweredMu[i - j - 1] * pow(fabs(presi[i - j - 1] - v) + c * (presi[i - j - 1] - v), d[1]); //adds the p-part
      for(j = 0; j < q; j++) poweredMu[i] += b[j] * poweredMu[i - j - 1]; //adds the q-part
      
      pmu[i] = pow(poweredMu[i], 1/d[0]);
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpq, N);
  } while (stopIndex != N);
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_AACDcall----------------------------//

//START---getLL_AACDcallEx----------------------------//
SEXP getLL_AACDcallEx(
    SEXP x,
    SEXP z, //the exogenous regressors
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], q = INTEGER(order)[1];
  int maxpq = max(p, q);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpq;
  
  double *px, *pz;
  int *pnewDay;
  px = REAL(x);
  pz = REAL(z);
  pnewDay = INTEGER(newDay);
  
  int N = length(x), NnewDays = length(newDay);
  int k = (int)(length(z)/N); //the number of exogenous variables 
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  
  double poweredMu[N]; //the mus to the d[0] power
  double c = REAL(par)[1 + 2 * p + q];
  double v = REAL(par)[2 + 2 * p + q];
  double d[2] = {REAL(par)[p + q + 3], REAL(par)[p + q + 4]};
  double a[p], b[q];
  
  for(j = 0; j < p; j++) a[j] = REAL(par)[j + 1];
  for(j = 0; j < q; j++) b[j] = REAL(par)[j + 1 + p];
  
  do {
    //in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
    for (i = startIndex; i < stopIndex; i++){
      pmu[i] = REAL(mean)[0];
      presi[i] = px[i]/pmu[i];
      poweredMu[i] = pow(REAL(mean)[0], d[0]);
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    for(i = startIndex; i < stopIndex; i++){
      poweredMu[i] = REAL(par)[0]; //adds the constant
      for(j = 0; j < p; j++) poweredMu[i] += a[j] * poweredMu[i - j - 1] * pow(fabs(presi[i - j - 1] - v) + c * (presi[i - j - 1] - v), d[1]); //adds the p-part
      for(j = 0; j < q; j++) poweredMu[i] += b[j] * poweredMu[i - j - 1]; //adds the q-part
      for(j = 0; j <= k - 1; j++) poweredMu[i] += REAL(par)[j + 1 + q + p] * pz[i + j * N]; //adds the exogenous variables			
      
      pmu[i] = pow(poweredMu[i], 1/d[0]);
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpq, N);
  } while (stopIndex != N);
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_AACDcallEx----------------------------//

//START---getLL_TACDcall----------------------------//
SEXP getLL_TACDcall(
    SEXP x,
    SEXP threshVar,
    SEXP threshType, // if 0, 'threshVar' is used as a threshold variable, if 1, the lagged cond. dur. mu is used
    SEXP bp,
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], q = INTEGER(order)[1];
  int maxpq = max(p, q);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpq;
  
  double *px;
  int *pnewDay;
  px = REAL(x);
  pnewDay = INTEGER(newDay);
  
  
  double *pthreshVar, tempThreshVar = 0, *pbp;
  pthreshVar = REAL(threshVar);
  pbp = REAL(bp);
  int currentThreshold = 0, J = length(bp) + 1;
  
  int N = length(x), NnewDays = length(newDay);
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  
  do {
    //at the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
    for(i = startIndex; i < stopIndex; i++){
      pmu[i] = REAL(mean)[0];
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    
    for(i = startIndex; i < stopIndex; i++){
      
      // gets the current regime:
      if(INTEGER(threshType)[0] == 0) tempThreshVar = pthreshVar[i - 1];
      if(INTEGER(threshType)[0] == 1) tempThreshVar = pmu[i - 1];
      currentThreshold = 0;
      if(tempThreshVar > pbp[0]){
        for(int i1 = 1; i1 < length(bp); i1++){
          if(tempThreshVar <= pbp[i1]){
            currentThreshold = i1;
            break;
          }
        }
        if(tempThreshVar > pbp[length(bp) - 1]) currentThreshold = length(bp);
      }
            
      pmu[i] = REAL(par)[0 + currentThreshold]; //adds the constant
      for(j = 1; j <= p; j++) pmu[i] += REAL(par)[J + currentThreshold * p + j - 1] * px[i - j]; //adds the lagged durations
      for(j = 1; j <= q; j++) pmu[i] += REAL(par)[(p + 1) * J + currentThreshold * q + j - 1] * pmu[i - j]; //adds the lagged mus
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpq, N);
  } while (stopIndex != N);
  
  
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_TACDcall----------------------------//

//START---getLL_TACDcallEx----------------------------//
SEXP getLL_TACDcallEx(
    SEXP x,
    SEXP z, //the exogenous regressors
    SEXP threshVar,
    SEXP threshType, // if 0, 'threshVar' is used as a threshold variable, if 1, the lagged cond. dur. mu is used
    SEXP bp,
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], q = INTEGER(order)[1];
  int maxpq = max(p, q);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpq;
  
  double *px, *pz;
  int *pnewDay;
  px = REAL(x);
  pz = REAL(z);
  pnewDay = INTEGER(newDay);
  
  
  double *pthreshVar, tempThreshVar = 0, *pbp;
  pthreshVar = REAL(threshVar);
  pbp = REAL(bp);
  int currentThreshold = 0, J = length(bp) + 1;
  
  int N = length(x), NnewDays = length(newDay);
  int k = (int)(length(z)/N); //the number of exogenous variables 
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  
  do {
    //at the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
    for(i = startIndex; i < stopIndex; i++){
      pmu[i] = REAL(mean)[0];
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    
    for(i = startIndex; i < stopIndex; i++){
      
      // gets the current regime:
      if(INTEGER(threshType)[0] == 0) tempThreshVar = pthreshVar[i - 1];
      if(INTEGER(threshType)[0] == 1) tempThreshVar = pmu[i - 1];
      currentThreshold = 0;
      if(tempThreshVar > pbp[0]){
        for(int i1 = 1; i1 < length(bp); i1++){
          if(tempThreshVar <= pbp[i1]){
            currentThreshold = i1;
            break;
          }
        }
        if(tempThreshVar > pbp[length(bp) - 1]) currentThreshold = length(bp);
      }
      
      pmu[i] = REAL(par)[currentThreshold]; //adds the constant
      for(j = 1; j <= p; j++) pmu[i] += REAL(par)[J + currentThreshold * p + j - 1] * px[i - j]; //adds the lagged durations
      for(j = 1; j <= q; j++) pmu[i] += REAL(par)[(p + 1) * J + currentThreshold * q + j - 1] * pmu[i - j]; //adds the lagged mus
      
      for(j = 0; j <= k - 1; j++) pmu[i] += REAL(par)[length(par) - k - 1 + j] * pz[i + j * N]; //adds the exogenous variables
      
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpq, N);
  } while (stopIndex != N);
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_TACDcallEx----------------------------//

//START---getLL_TAMACDcall----------------------------//
SEXP getLL_TAMACDcall(
    SEXP x,
    SEXP threshVar,
    SEXP threshType, // if 0, 'threshVar' is used as a threshold variable, if 1, the lagged cond. dur. mu is used
    SEXP bp,
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], r = INTEGER(order)[1], q = INTEGER(order)[2];
  int maxpqr = max(p, q); maxpqr = max(maxpqr, r);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpqr;
  
  double *px;
  int *pnewDay;
  px = REAL(x);
  pnewDay = INTEGER(newDay);
  
  
  double *pthreshVar, tempThreshVar = 0, *pbp;
  pthreshVar = REAL(threshVar);
  pbp = REAL(bp);
  int currentThreshold = 0, J = length(bp) + 1;
  
  int N = length(x), NnewDays = length(newDay);
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  
  do {
    //at the start of the sample or at the start of a new day, the maxpqr mus are set to the mean:
    for(i = startIndex; i < stopIndex; i++){
      pmu[i] = REAL(mean)[0];
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    
    for(i = startIndex; i < stopIndex; i++){
      
      // gets the current regime:
      if(INTEGER(threshType)[0] == 0) tempThreshVar = pthreshVar[i - 1];
      if(INTEGER(threshType)[0] == 1) tempThreshVar = pmu[i - 1];
      currentThreshold = 0;
      if(tempThreshVar > pbp[0]){
        for(int i1 = 1; i1 < length(bp); i1++){
          if(tempThreshVar <= pbp[i1]){
            currentThreshold = i1;
            break;
          }
        }
        if(tempThreshVar > pbp[length(bp) - 1]) currentThreshold = length(bp);
      }	  
	  
      pmu[i] = REAL(par)[0 + currentThreshold]; //adds the constant
      for(j = 1; j <= p; j++) pmu[i] += REAL(par)[J + currentThreshold * p + j - 1] * px[i - j]; //adds the lagged durations
        	  
	  for(j = 1; j <= r; j++) pmu[i] += REAL(par)[(p + 1) * J + currentThreshold * r + j - 1] * presi[i - j]; //adds the lagged residuals	   
	  
	  for(j = 1; j <= q; j++) pmu[i] += REAL(par)[(p + r + 1) * J + currentThreshold * q + j - 1] * pmu[i - j]; //adds the lagged mus
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpqr, N);
  } while (stopIndex != N);
  
  
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_TAMACDcall----------------------------//

//START---getLL_TAMACDcallEx----------------------------//
SEXP getLL_TAMACDcallEx(
    SEXP x,
    SEXP z, //the exogenous regressors
    SEXP threshVar,
    SEXP threshType, // if 0, 'threshVar' is used as a threshold variable, if 1, the lagged cond. dur. mu is used
    SEXP bp,
    SEXP par,
    SEXP order,
    SEXP mean,
    SEXP dist,
    SEXP distPara,
    SEXP newDay,
    SEXP forceErrExpec){
  
  
  int p = INTEGER(order)[0], r = INTEGER(order)[1], q = INTEGER(order)[2];
  int maxpqr = max(p, q); maxpqr = max(maxpqr, r);
  int i = 0;
  int j = 1;
  int nextND=0;
  int startIndex = 0, stopIndex = maxpqr;
  
  double *px, *pz;
  int *pnewDay;
  px = REAL(x);
  pz = REAL(z);
  pnewDay = INTEGER(newDay);
  
  
  double *pthreshVar, tempThreshVar = 0, *pbp;
  pthreshVar = REAL(threshVar);
  pbp = REAL(bp);
  int currentThreshold = 0, J = length(bp) + 1;
  
  int N = length(x), NnewDays = length(newDay);
  int k = (int)(length(z)/N); //the number of exogenous variables 
  if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
  SEXP mu, resi;
  PROTECT(mu = NEW_NUMERIC(N));
  PROTECT(resi = NEW_NUMERIC(N));
  double *pmu, *presi;
  pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
  
  do {
    //at the start of the sample or at the start of a new day, the maxpqr mus are set to the mean:
    for(i = startIndex; i < stopIndex; i++){
      pmu[i] = REAL(mean)[0];
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
    else stopIndex = N;
    
    for(i = startIndex; i < stopIndex; i++){
      
      // gets the current regime:
      if(INTEGER(threshType)[0] == 0) tempThreshVar = pthreshVar[i - 1];
      if(INTEGER(threshType)[0] == 1) tempThreshVar = pmu[i - 1];
      currentThreshold = 0;
      if(tempThreshVar > pbp[0]){
        for(int i1 = 1; i1 < length(bp); i1++){
          if(tempThreshVar <= pbp[i1]){
            currentThreshold = i1;
            break;
          }
        }
        if(tempThreshVar > pbp[length(bp) - 1]) currentThreshold = length(bp);
      }	  
	  
      pmu[i] = REAL(par)[0 + currentThreshold]; //adds the constant
      for(j = 1; j <= p; j++) pmu[i] += REAL(par)[J + currentThreshold * p + j - 1] * px[i - j]; //adds the lagged durations
        	  
	  for(j = 1; j <= r; j++) pmu[i] += REAL(par)[(p + 1) * J + currentThreshold * r + j - 1] * presi[i - j]; //adds the lagged residuals	   
	  
	  for(j = 1; j <= q; j++) pmu[i] += REAL(par)[(p + r + 1) * J + currentThreshold * q + j - 1] * pmu[i - j]; //adds the lagged mus
	  
	  for(j = 0; j <= k - 1; j++) pmu[i] += REAL(par)[length(par) - k - 1 + j] * pz[i + j * N]; //adds the exogenous variables
	  
      presi[i] = px[i]/pmu[i];
    }
    startIndex = stopIndex;
    stopIndex = min(stopIndex + maxpqr, N);
  } while (stopIndex != N);
  
  
  
  SEXP list, LL;
  PROTECT(LL = NEW_NUMERIC(1));
  PROTECT(list = NEW_LIST(3));
  
  SET_VECTOR_ELT(list, 0, mu);
  SET_VECTOR_ELT(list, 1, resi);
  
  REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
  SET_VECTOR_ELT(list, 2, LL);
  
  UNPROTECT(4);
  return list;
}
//END---getLL_TAMACDcallEx----------------------------//


//START---getLL_AMACDcall----------------------------//
SEXP getLL_AMACDcall(
		SEXP x,
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){


	int p = INTEGER(order)[0], r = INTEGER(order)[1], q = INTEGER(order)[2];
	int maxpqr = max(max(p, q), r);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpqr;

	double *px;
	int *pnewDay;
	px = REAL(x);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(mean)[0];
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(par)[0]; //adds the constant
			for(j = 1; j <= p; j++) pmu[i] += REAL(par)[j] * px[i - j]; //adds the lagged durations
			for(j = 1; j <= r; j++) pmu[i] += REAL(par)[j + p] * presi[i - 1];  //adds the lagged residuals
			for(j = 1; j <= q; j++) pmu[i] += REAL(par)[j + r + p] * pmu[i - j]; //adds the lagged mus
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpqr, N);
	} while (stopIndex != N);	
	
	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_AMACDcall----------------------------//

//START---getLL_AMACDcallEx----------------------------//
SEXP getLL_AMACDcallEx(
		SEXP x,
		SEXP z, //the exogenous regressors
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){


	int p = INTEGER(order)[0], r = INTEGER(order)[1], q = INTEGER(order)[2];
	int maxpqr = max(max(p, q), r);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpqr;

	double *px, *pz;
	int *pnewDay;
	px = REAL(x);
	pz = REAL(z);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	int k = (int)(length(z)/N); //the number of exogenous variables 
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(mean)[0];
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			pmu[i] = REAL(par)[0]; //adds the constant
			for(j = 1; j <= p; j++) pmu[i] += REAL(par)[j] * px[i - j]; //adds the lagged durations
			for(j = 1; j <= r; j++) pmu[i] += REAL(par)[j + p] * presi[i - 1];  //adds the lagged residuals
			for(j = 1; j <= q; j++) pmu[i] += REAL(par)[j + r + p] * pmu[i - j]; //adds the lagged mus
			for(j = 0; j <= k - 1; j++) pmu[i] += REAL(par)[j + 1 + q + p] * pz[i + j * N]; //adds the exogenous variables			
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpqr, N);
	} while (stopIndex != N);	
	
	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_AMACDcallEx----------------------------//



//START---getLL_LACD2call----------------------------//
SEXP getLL_LACD2call(
		SEXP x,
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){


	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	int maxpq = max(p, q);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpq;

	double *px;
	int *pnewDay;
	px = REAL(x);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
	double logMu[N];

	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			logMu[i] = log(REAL(mean)[0]);
			pmu[i] = exp(logMu[i]);
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			logMu[i] = REAL(par)[0]; //adds the constant
			for(j = 1; j <= p; j++) logMu[i] += REAL(par)[j] * presi[i - j]; //adds the lagged durations
			for(j = 1; j <= q; j++) logMu[i] += REAL(par)[j + p] * logMu[i - j]; //adds the lagged mus
			pmu[i] = exp(logMu[i]);
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, N);
	} while (stopIndex != N);
	
	
	
	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_LACD2call----------------------------//

//START---getLL_LACD2callEx----------------------------//
SEXP getLL_LACD2callEx(
		SEXP x,
		SEXP z, //the exogenous regressors
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP forceErrExpec){


	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	int maxpq = max(p, q);
	int i = 0;
	int j = 1;
	int nextND=0;
	int startIndex = 0, stopIndex = maxpq;

	double *px, *pz;
	int *pnewDay;
	px = REAL(x);
	pz = REAL(z);
	pnewDay = INTEGER(newDay);

	int N = length(x), NnewDays = length(newDay);
	int k = (int)(length(z)/N); //the number of exogenous variables 
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
	SEXP mu, resi;
	PROTECT(mu = NEW_NUMERIC(N));
	PROTECT(resi = NEW_NUMERIC(N));
	double *pmu, *presi;
	pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);
	double logMu[N];

	do {
		//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			logMu[i] = log(REAL(mean)[0]);
			pmu[i] = exp(logMu[i]);
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;
		for(i = startIndex; i < stopIndex; i++){
			logMu[i] = REAL(par)[0]; //adds the constant
			for(j = 1; j <= p; j++) logMu[i] += REAL(par)[j] * presi[i - j]; //adds the lagged durations
			for(j = 1; j <= q; j++) logMu[i] += REAL(par)[j + p] * logMu[i - j]; //adds the lagged mus
			for(j = 0; j <= k - 1; j++) logMu[i] += REAL(par)[j + 1 + q + p] * pz[i + j * N]; //adds the exogenous variables			
			pmu[i] = exp(logMu[i]);
			presi[i] = px[i]/pmu[i];
		}
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, N);
	} while (stopIndex != N);
	
	
	
	SEXP list, LL;
	PROTECT(LL = NEW_NUMERIC(1));
	PROTECT(list = NEW_LIST(3));

	SET_VECTOR_ELT(list, 0, mu);
	SET_VECTOR_ELT(list, 1, resi);

	REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));
	SET_VECTOR_ELT(list, 2, LL);

	UNPROTECT(4);
	return list;
}
//END---getLL_LACD2callEx----------------------------//

//START---getLL_SNIACDcall----------------------------//
SEXP getLL_SNIACDcall(		SEXP x,
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP bp,
		SEXP forceErrExpec){

int p = INTEGER(order)[0], q = INTEGER(order)[1];
int maxpq = max(p, q);
int i = 0;
int j = 1;
int k = 0;
int nextND=0;
int startIndex = 0, stopIndex = maxpq;

double *px;
int *pnewDay;
px = REAL(x);
pnewDay = INTEGER(newDay);

int N = length(x), NnewDays = length(newDay);
if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
SEXP mu, resi;
PROTECT(mu = NEW_NUMERIC(N));
PROTECT(resi = NEW_NUMERIC(N));
double *pmu, *presi;
pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

int M = length(bp);
double a[p], b[q], c[M+1]; //parameters

//puts the parameter pointers to local variables for easier interpretation:
for(j = 0; j < (M+1); j++){
	c[j] = REAL(par)[j+1];
}
a[0] = 0;
for(j = 1; j < (p-1); j++){
	a[j] = REAL(par)[j+1+M];
}
for(j = 0; j < q; j++){
	b[j] = REAL(par)[j+1+M+p];
}



do {
	//at the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
	for (i = startIndex; i < stopIndex; i++){
		pmu[i] = REAL(mean)[0];
		presi[i] = px[i]/pmu[i];
	}
	startIndex = stopIndex;
	if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
	else stopIndex = N;
	for(i = startIndex; i < stopIndex; i++){
		pmu[i] = REAL(par)[0]; //adds the constant

		for(j = 0; j < p; j++){
			pmu[i] += (c[0] + a[j]) * presi[i-j-1];
			k = 0;
			while(presi[i-j-1] >= REAL(bp)[k] && k < M){
				pmu[i] += (c[k+1] + a[j])*(presi[i-j-1] - REAL(bp)[k]);
				k++;
			}
		}

		for(j = 0; j < q; j++) pmu[i] += b[j]* pmu[i - j - 1]; //adds the q-part
		presi[i] = px[i]/pmu[i]; //computes the residual
	}
	startIndex = stopIndex;
	stopIndex = min(stopIndex + maxpq, N);
} while (stopIndex != N);

SEXP list, LL;
PROTECT(LL = NEW_NUMERIC(1));
PROTECT(list = NEW_LIST(3));

REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));

SET_VECTOR_ELT(list, 0, mu);
SET_VECTOR_ELT(list, 1, resi);
SET_VECTOR_ELT(list, 2, LL);

UNPROTECT(4);
return list;
}
//END---getLL_SNIACDcall----------------------------//

//START---getLL_SNIACDcallEx----------------------------//
SEXP getLL_SNIACDcallEx(
		SEXP x,
		SEXP z, //the exogenous regressors
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP bp,
		SEXP forceErrExpec){

int p = INTEGER(order)[0], q = INTEGER(order)[1];
int maxpq = max(p, q);
int i = 0;
int j = 1;
int k = 0;
int nextND=0;
int startIndex = 0, stopIndex = maxpq;

	double *px, *pz;
int *pnewDay;
px = REAL(x);
pz = REAL(z);
pnewDay = INTEGER(newDay);

int N = length(x), NnewDays = length(newDay);
int g = (int)(length(z)/N); //the number of exogenous variables 
if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
SEXP mu, resi;
PROTECT(mu = NEW_NUMERIC(N));
PROTECT(resi = NEW_NUMERIC(N));
double *pmu, *presi;
pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

int M = length(bp);
double a[p], b[q], c[M+1]; //parameters

//puts the parameter pointers to local variables for easier interpretation:
for(j = 0; j < (M+1); j++){
	c[j] = REAL(par)[j+1];
}
a[0] = 0;
for(j = 1; j < (p-1); j++){
	a[j] = REAL(par)[j+1+M];
}
for(j = 0; j < q; j++){
	b[j] = REAL(par)[j+1+M+p];
}



do {
	//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
	for (i = startIndex; i < stopIndex; i++){
		pmu[i] = REAL(mean)[0];
		presi[i] = px[i]/pmu[i];
	}
	startIndex = stopIndex;
	if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
	else stopIndex = N;
	for(i = startIndex; i < stopIndex; i++){
		pmu[i] = REAL(par)[0]; //adds the constant

		for(j = 0; j < p; j++){
			pmu[i] += (c[0] + a[j]) * presi[i-j-1];
			k = 0;
			while(presi[i-j-1] >= REAL(bp)[k] && k < M){
				pmu[i] += (c[k+1] + a[j])*(presi[i-j-1] - REAL(bp)[k]);
				k++;
			}
		}

		for(j = 0; j < q; j++) pmu[i] += b[j]* pmu[i - j - 1]; //adds the q-part
		for(j = 0; j <= g - 1; j++) pmu[i] += REAL(par)[j + 1 + q + p] * pz[i + j * N]; //adds the exogenous variables			
			
		presi[i] = px[i]/pmu[i]; //computes the residual
	}
	startIndex = stopIndex;
	stopIndex = min(stopIndex + maxpq, N);
} while (stopIndex != N);

SEXP list, LL;
PROTECT(LL = NEW_NUMERIC(1));
PROTECT(list = NEW_LIST(3));

REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), 
INTEGER(forceErrExpec));

SET_VECTOR_ELT(list, 0, mu);
SET_VECTOR_ELT(list, 1, resi);
SET_VECTOR_ELT(list, 2, LL);

UNPROTECT(4);
return list;
}
//END---getLL_SNIACDcallEx----------------------------//

//START---getLL_logSNIACDcall----------------------------//
SEXP getLL_logSNIACDcall(		SEXP x,
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP bp,
		SEXP forceErrExpec){

int p = INTEGER(order)[0], q = INTEGER(order)[1];
int maxpq = max(p, q);
int i = 0;
int j = 1;
int k = 0;
int nextND=0;
int startIndex = 0, stopIndex = maxpq;

double *px;
int *pnewDay;
px = REAL(x);
pnewDay = INTEGER(newDay);

int N = length(x), NnewDays = length(newDay);
if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
SEXP mu, resi;
PROTECT(mu = NEW_NUMERIC(N));
PROTECT(resi = NEW_NUMERIC(N));
double logMu[N];
double *pmu, *presi;
pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

int M = length(bp);
double a[p], b[q], c[M+1]; //parameters

//puts the parameter pointers to local variables for easier interpretation:
for(j = 0; j < (M+1); j++){
	c[j] = REAL(par)[j+1];
}
a[0] = 0;
for(j = 1; j < (p-1); j++){
	a[j] = REAL(par)[j+1+M];
}
for(j = 0; j < q; j++){
	b[j] = REAL(par)[j+1+M+p];
}

do {
	//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
	for (i = startIndex; i < stopIndex; i++){
		
		pmu[i] = REAL(mean)[0];
		logMu[i] = log(pmu[i]);
		presi[i] = px[i]/pmu[i];
	}
	startIndex = stopIndex;
	if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
	else stopIndex = N;
	for(i = startIndex; i < stopIndex; i++){
		logMu[i] = REAL(par)[0]; //adds the constant

		for(j = 0; j < p; j++){
			logMu[i] += (c[0] + a[j]) * presi[i-j-1];
			k = 0;
			while(presi[i-j-1] >= REAL(bp)[k] && k < M){
				logMu[i] += (c[k+1] + a[j])*(presi[i-j-1] - REAL(bp)[k]);
				k++;
			}
		}

		for(j = 0; j < q; j++) logMu[i] += b[j]* logMu[i - j - 1]; //adds the q-part
		
		pmu[i] = exp(logMu[i]);
		presi[i] = px[i]/pmu[i]; //computes the residual
	}
	startIndex = stopIndex;
	stopIndex = min(stopIndex + maxpq, N);
} while (stopIndex != N);

SEXP list, LL;
PROTECT(LL = NEW_NUMERIC(1));
PROTECT(list = NEW_LIST(3));

REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), INTEGER(forceErrExpec));

SET_VECTOR_ELT(list, 0, mu);
SET_VECTOR_ELT(list, 1, resi);
SET_VECTOR_ELT(list, 2, LL);

UNPROTECT(4);
return list;
}
//END---getLL_logSNIACDcall----------------------------//


//START---getLL_logSNIACDcallEx----------------------------//
SEXP getLL_logSNIACDcallEx(
		SEXP x,
		SEXP z, //the exogenous regressors
		SEXP par,
		SEXP order,
		SEXP mean,
		SEXP dist,
		SEXP distPara,
		SEXP newDay,
		SEXP bp,
		SEXP forceErrExpec){

int p = INTEGER(order)[0], q = INTEGER(order)[1];
int maxpq = max(p, q);
int i = 0;
int j = 1;
int k = 0;
int nextND=0;
int startIndex = 0, stopIndex = maxpq;

double *px, *pz;
int *pnewDay;
px = REAL(x);
pz = REAL(z);
pnewDay = INTEGER(newDay);

int N = length(x), NnewDays = length(newDay);
int g = (int)(length(z)/N); //the number of exogenous variables 
if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0;
SEXP mu, resi;
PROTECT(mu = NEW_NUMERIC(N));
PROTECT(resi = NEW_NUMERIC(N));
double logMu[N];
double *pmu, *presi;
pmu = NUMERIC_POINTER(mu); presi = NUMERIC_POINTER(resi);

int M = length(bp);
double a[p], b[q], c[M+1]; //parameters

//puts the parameter pointers to local variables for easier interpretation:
for(j = 0; j < (M+1); j++){
	c[j] = REAL(par)[j+1];
}
a[0] = 0;
for(j = 1; j < (p-1); j++){
	a[j] = REAL(par)[j+1+M];
}
for(j = 0; j < q; j++){
	b[j] = REAL(par)[j+1+M+p];
}

do {
	//in the start of the sample or at the start of a new day, the maxpq mus are set to the mean:
	for (i = startIndex; i < stopIndex; i++){
		
		pmu[i] = REAL(mean)[0];
		logMu[i] = log(pmu[i]);
		presi[i] = px[i]/pmu[i];
	}
	startIndex = stopIndex;
	if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
	else stopIndex = N;
	for(i = startIndex; i < stopIndex; i++){
		logMu[i] = REAL(par)[0]; //adds the constant

		for(j = 0; j < p; j++){
			logMu[i] += (c[0] + a[j]) * presi[i-j-1];
			k = 0;
			while(presi[i-j-1] >= REAL(bp)[k] && k < M){
				logMu[i] += (c[k+1] + a[j])*(presi[i-j-1] - REAL(bp)[k]);
				k++;
			}
		}

		for(j = 0; j < q; j++) logMu[i] += b[j]* logMu[i - j - 1]; //adds the q-part
		for(j = 0; j <= g - 1; j++) logMu[i] += REAL(par)[j + 1 + q + p] * pz[i + j * N]; //adds the exogenous variables			
			
		pmu[i] = exp(logMu[i]);		
		presi[i] = px[i]/pmu[i]; //computes the residual
	}
	startIndex = stopIndex;
	stopIndex = min(stopIndex + maxpq, N);
} while (stopIndex != N);

SEXP list, LL;
PROTECT(LL = NEW_NUMERIC(1));
PROTECT(list = NEW_LIST(3));

REAL(LL)[0] = getLL_dist(px, pmu, presi, &N, INTEGER(dist), REAL(distPara), 
INTEGER(forceErrExpec));

SET_VECTOR_ELT(list, 0, mu);
SET_VECTOR_ELT(list, 1, resi);
SET_VECTOR_ELT(list, 2, LL);

UNPROTECT(4);
return list;
}
//END---getLL_logSNIACDcallEx----------------------------//
