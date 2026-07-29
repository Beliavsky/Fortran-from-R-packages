#include "header.h" 
#include <Rinternals.h>
#include <Rdefines.h>


//START---sim_BACD----------------------------//
SEXP sim_BACD(SEXP N,
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
  double poweredMu[INTEGER(N)[0] + INTEGER(Nburn)[0]];
  double xTemp[INTEGER(N)[0] + INTEGER(Nburn)[0]];
  double d[2] = {xpar[1 + xorder[0] + xorder[1]], 
                 xpar[2 + xorder[0] + xorder[1]]};
  
  for (i = 0; i < Nstart; i++){
    xTemp[i] = xstartX[i];
    mu[i] = xstartMu[i];
    poweredMu[i] = pow(mu[i], d[0]);
  }
  for(i = Nstart; i < INTEGER(N)[0] + INTEGER(Nburn)[0]; i++){
    
      poweredMu[i] = xpar[0]; //adds the constant
      for(j = 1; j <= xorder[0]; j++) poweredMu[i] += xpar[j] * pow(xe[i-j], d[1]); //adds the lagged durations
      for(j = 1; j <= xorder[1]; j++) poweredMu[i] += xpar[j+xorder[0]] * poweredMu[i-j]; //adds the lagged mus
      
      xTemp[i] = pow(poweredMu[i], 1 / d[0]) * xe[i];
    
  }
  
  PROTECT(dur = NEW_NUMERIC(INTEGER(N)[0]));
  xdur = NUMERIC_POINTER (dur) ;
  
  for(i =  0; i < INTEGER(N)[0]; i++){
    xdur[i] = xTemp[i + INTEGER(Nburn)[0]];
  }
  
  UNPROTECT(7);
  return(dur);
}
//END---sim_BACD----------------------------//
