#include "header.h" 
#include <Rinternals.h>
#include <Rdefines.h>

//START---sim_TACD----------------------------//
SEXP sim_TACD(SEXP N,
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
  int regime[INTEGER(N)[0] + INTEGER(Nburn)[0]];
  
  int p = INTEGER(order)[0], q = INTEGER(order)[1];
  
  int M = length(bp); // M is the number of breakpoints
  double omega[M + 1], alpha[p * (M + 1)], beta[q * (M + 1)]; //parameters
  
  for(j = 0; j < (M + 1); j++) omega[j] = xpar[j];
  for(j = 0; j < (M + 1) * p; j++) alpha[j] = xpar[(M + 1) + j];
  for(j = 0; j < (M + 1) * q; j++) beta[j] = xpar[(p + 1) * (M + 1) + j];
  
  for (i = 0; i < Nstart; i++){
    mu[i] = xstartMu[i];
    xTemp[i] = xstartX[i];
  }
  for(i = Nstart; i < INTEGER(N)[0] + INTEGER(Nburn)[0]; i++){
    
    // gets the current regime:
    regime[i] = 0;
    double tempThreshVar = xTemp[i - 1]; 
    if(tempThreshVar >= REAL(bp)[0]){
      for(int bpIndex = 1; bpIndex < M; bpIndex++){
        if(tempThreshVar < REAL(bp)[bpIndex]){
          regime[i] = bpIndex;
          break;
        }
      }
      if(tempThreshVar >= REAL(bp)[M - 1]) regime[i] = M;
    }
    
    mu[i] = omega[regime[i]]; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < p; j++) mu[i] += alpha[j * (M + 1) + regime[i]] * xTemp[i - j - 1];
    
    //adds the "q"-part:
    for(int j = 0; j < q; j++) mu[i] += beta[j * (M + 1) + regime[i]] * mu[i - j - 1]; 
    
    
    xTemp[i] = mu[i] * xe[i];
  }
  
  PROTECT(dur = NEW_NUMERIC(INTEGER(N)[0])); Nprotect++;
  xdur = NUMERIC_POINTER (dur);
  
  // puts the simulated durations, less the first Nburn durations, in 'dur':
  for(i =  0; i < INTEGER(N)[0]; i++){
    xdur[i] = xTemp[i + INTEGER(Nburn)[0]];
  }
  
  UNPROTECT(Nprotect);
  return(dur);
}
//END---sim_TACD----------------------------//
