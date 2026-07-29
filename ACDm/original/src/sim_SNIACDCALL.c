#include "header.h" 
#include <Rinternals.h>
#include <Rdefines.h>

//START---sim_SNIACD----------------------------//
SEXP sim_SNIACD(SEXP N,
		SEXP par,
		SEXP order,
		SEXP startX,
		SEXP startMu,
		SEXP e,
		SEXP Nburn,
		SEXP bp){

	int i = 0;
	int j = 1;
	SEXP dur;
	int Nprotect = 0;
	PROTECT(N = AS_INTEGER(N)); Nprotect++;
	double *xpar, *xstartX, *xstartMu, *xe, *xdur;
	PROTECT(par = AS_NUMERIC(par)); Nprotect++;
	PROTECT(order = AS_INTEGER(order)); Nprotect++;
	PROTECT(startX = AS_NUMERIC(startX)); Nprotect++;
	PROTECT(e = AS_NUMERIC(e)); Nprotect++;
	PROTECT(Nburn = AS_INTEGER(Nburn)); Nprotect++;
	PROTECT(bp = AS_NUMERIC(bp)); Nprotect++;
	
	int Nstart = LENGTH(startX);
	xpar = NUMERIC_POINTER(par); 
	xstartX = NUMERIC_POINTER(startX); xstartMu = NUMERIC_POINTER(startMu);
	xe = NUMERIC_POINTER(e);
	
	double mu[INTEGER(N)[0] + INTEGER(Nburn)[0]];
	double xTemp[INTEGER(N)[0] + INTEGER(Nburn)[0]];
	
	int p = INTEGER(order)[0], q = INTEGER(order)[1];
	
	int M = length(bp);
	double a[p], b[q], c[M+1]; //parameters
	
	//puts the parameter pointers to local variables for easier interpretation:
	for(j = 0; j < (M+1); j++){
	  c[j] = xpar[j+1];
	}
	a[0] = 0;
	for(j = 1; j < (p-1); j++){
	  a[j] = xpar[j+1+M];
	}
	for(j = 0; j < q; j++){
	  b[j] = xpar[j+1+M+p];
	}

	for (i = 0; i < Nstart; i++){
		mu[i] = xstartMu[i];
		xTemp[i] = xstartX[i];
	}
	for(i = Nstart; i < INTEGER(N)[0] + INTEGER(Nburn)[0]; i++){
		mu[i] = xpar[0]; //adds the constant

	  //adds the spline part:
		for(j = 0; j < p; j++){
		  double holder = 0;
		  holder += c[0] * xe[i-j-1];
		  
		  for(int k = 1; k < M + 1; k++){
		    if(xe[i - j - 1] > REAL(bp)[k - 1]) holder += c[k] * (xe[i - j - 1] - REAL(bp)[k - 1]);
		    else break;
		  }
		  
		  if(j > 0) holder *= a[j - 1];
		  mu[i] += holder;
		}
		
		//adds the "q" part:
		for(j = 0; j < q; j++) mu[i] += b[j] * mu[i - j - 1];
		
		
		xTemp[i] = mu[i] * xe[i];
	}

	PROTECT(dur = NEW_NUMERIC(INTEGER(N)[0])); Nprotect++;
	xdur = NUMERIC_POINTER (dur) ;
	
	for(i =  0; i < INTEGER(N)[0]; i++){
		xdur[i] = xTemp[i + INTEGER(Nburn)[0]];
	}
	
	UNPROTECT(Nprotect);
	return(dur);
}
//END---sim_SNIACD----------------------------//
