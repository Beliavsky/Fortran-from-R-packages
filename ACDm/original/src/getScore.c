#include "header.h" 
#include <Rinternals.h>
#include <Rdefines.h>

//START---getScoreACDExp----------------------------//
//    calculates the expected score and hessian
//    returns the expected score and derivative of the mean 
//    for each observation, and the summed expected hessian,
//    as well as the outer product of the scores
SEXP getScoreACDExp(
		SEXP x,
		SEXP mu,
		SEXP par,
		SEXP order,
		SEXP newDay){

	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	int maxpq = max(p, q);
	int i = 0, j = 0, v = 0, N = length(x), row = 0, col;
	int Npara = INTEGER(order)[0] + INTEGER(order)[1] + 1;
	int startIndex = 0, stopIndex = maxpq, nextND=0, NnewDays = length(newDay);
	int *pnewDay; pnewDay = INTEGER(newDay);
	if(NnewDays == 1 && pnewDay[0] == 0) NnewDays = 0; //in case there are no new days

	SEXP dmydtheta;
	PROTECT(dmydtheta = allocMatrix(REALSXP, N, Npara));
	double *dmydthetaptr; dmydthetaptr = REAL(dmydtheta);

	SEXP dLdtheta;
	PROTECT(dLdtheta = allocMatrix(REALSXP, N, Npara));
	double *dLdthetaptr; dLdthetaptr = REAL(dLdtheta);

	SEXP hessian;
	PROTECT(hessian = allocMatrix(REALSXP, Npara, Npara));
	double *hessianptr; hessianptr = REAL(hessian);

	SEXP OPscore; //sum of the outer product of the scores
	PROTECT(OPscore = allocMatrix(REALSXP, Npara, Npara));
	double *OPscoreptr; OPscoreptr = REAL(OPscore);

	double *px; px = REAL(x);
	double *pmu; pmu = REAL(mu);


	//initiates the hessian and OPscore:
	for(j = 0; j < Npara * Npara; j++){
		hessianptr[j] = 0;
		OPscoreptr[j] = 0;
	}


	do{
		//fills the first Npara number of dmydtheta with 0's, as the log likelihood 
		//for these observations are unchanged by the parameters:
		for(i = startIndex; i < stopIndex; i++) {
			for(j = 0; j < Npara; j++){
				dmydthetaptr[i + j * N] = 0;
				dLdthetaptr[i + j * N] = 0;
			}
		}
		
		startIndex = stopIndex;
		if(nextND < NnewDays) stopIndex = pnewDay[nextND++] - 1;
		else stopIndex = N;

		//fills the rest of dmydtheta:
		for (i = startIndex; i < stopIndex; i++) {

			//for omega (the constant)
			dmydthetaptr[i] = 1;
			for(j = 1; j <= q; j++) dmydthetaptr[i] +=
					REAL(par)[j + p] * dmydthetaptr[i - j];

			//for alpha
			for(v = 1; v <= p; v++){
				dmydthetaptr[i + v * N] = px[i - v];
				for(j = 1; j <= q; j++) dmydthetaptr[i + v * N] +=
						REAL(par)[j + p] * dmydthetaptr[i - j + v * N];
			}

			//for beta
			for(v = 1; v <= q; v++){
				dmydthetaptr[i + (v + p) * N] = pmu[i - v];
				for(j = 1; j <= q; j++) dmydthetaptr[i + (v + p) * N] +=
						REAL(par)[j + p] * dmydthetaptr[i - j + (v + p) * N];
			}

			//calculates the derivatives of the log likelihood:
			for(j = 0; j < Npara; j++){
				dLdthetaptr[i + j * N] =
						dmydthetaptr[i + j * N] * (-1 / pmu[i] + px[i] / (pmu[i] * pmu[i]) );
			}

			//calculates and adds the score outer product for observation i
			for(row = 0; row < Npara ; row++){ //row
				for(col = 0; col < Npara; col++){  //column
					OPscoreptr[row + col * Npara] += dLdthetaptr[i + row * N] * dLdthetaptr[i + col * N];
				}
			}

			//calculates and adds the hessian for observation i
			for(row = 0; row < Npara ; row++){ //row
				for(col = 0; col < Npara; col++){  //column
					hessianptr[row + col * Npara] -= pow(pmu[i], -2) * dmydthetaptr[i + row * N] * dmydthetaptr[i + col * N];
				}
			}
		}

		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, N);

	} while (stopIndex != N);

	SEXP list; PROTECT(list = NEW_LIST(4));

	SET_VECTOR_ELT(list, 0, dmydtheta);
	SET_VECTOR_ELT(list, 1, dLdtheta);
	SET_VECTOR_ELT(list, 2, hessian);
	SET_VECTOR_ELT(list, 3, OPscore);
	UNPROTECT(5);
	return list;
}
////END---getScoreACDExp----------------------------//
