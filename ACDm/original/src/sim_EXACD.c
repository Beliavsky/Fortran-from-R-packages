#include "header.h" 
#include <Rinternals.h>
#include <Rdefines.h>

//START---sim_EXACD----------------------------//
SEXP sim_EXACD(SEXP N,
               SEXP par,
               SEXP order,
               SEXP startX,
               SEXP startMu,
               SEXP e,
               SEXP Nburn){
  
  int i = 0;
  int j = 1;
  SEXP dur;
  PROTECT(N = AS_INTEGER(N));
  double *xpar, *xstartX, *xstartMu, *xe, *xdur;
  int *xorder;
  PROTECT(par = AS_NUMERIC(par));
  PROTECT(order = AS_INTEGER(order));
  PROTECT(startX = AS_NUMERIC(startX));
  PROTECT(e = AS_NUMERIC(e));
  PROTECT(Nburn = AS_INTEGER(Nburn));
  
  int Nstart = LENGTH(startX);
  xpar = NUMERIC_POINTER(par); xorder = INTEGER_POINTER(order);
  xstartX = NUMERIC_POINTER(startX); xstartMu = NUMERIC_POINTER(startMu);
  xe = NUMERIC_POINTER(e);
  
  double mu[INTEGER(N)[0] + INTEGER(Nburn)[0]];
  double LogMu[INTEGER(N)[0] + INTEGER(Nburn)[0]];
  double xTemp[INTEGER(N)[0] + INTEGER(Nburn)[0]];
  
  double a[xorder[0]], b[xorder[1]], d[xorder[0]];
  
  for(j = 0; j < xorder[0]; j++){
    a[j] = xpar[j + 1];
    d[j] = xpar[j + 1 + xorder[0]];
  }
  for(j = 0; j < xorder[1]; j++)
    b[j] = xpar[j + 1 + 2 * xorder[0]];
  
  // fills the first observations with their start values:
  for (i = 0; i < Nstart; i++){
    xTemp[i] = xstartX[i];
    mu[i] = xstartMu[i];
    LogMu[i] = log(mu[i]);
  }
  
  // computes the rest of the observations:
  for(i = Nstart; i < INTEGER(N)[0] + INTEGER(Nburn)[0]; i++){
    LogMu[i] = xpar[0]; //adds the constant
    //adds the lagged durations:
    for(j = 1; j <= xorder[0]; j++){ 
      LogMu[i] += a[j - 1] * xe[i - j]
                  + d[j - 1] * fabs(xe[i - j] - 1);
    }
    //adds the lagged mus:
    for(j = 1; j <= xorder[1]; j++){
      LogMu[i] += b[j - 1] * LogMu[i - j]; 
      xTemp[i] = exp(LogMu[i]) * xe[i];
    } 
  }
  
  PROTECT(dur = NEW_NUMERIC(INTEGER(N)[0]));
  xdur = NUMERIC_POINTER (dur) ;
  
  // fills 'dur' with only the last N observations:
  for(i =  0; i < INTEGER(N)[0]; i++){
    xdur[i] = xTemp[i + INTEGER(Nburn)[0]];
  }
  
  UNPROTECT(7);
  return(dur);
}
//END---sim_EXACD----------------------------//
