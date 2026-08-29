#include "header.h"

SEXP gibbs(SEXP n, SEXP np, SEXP thin, SEXP init,
	   SEXP propsd, SEXP f, SEXP env){

  int i,j,k,nr, l = 0;
  int nn = INTEGER(n)[0], nnp = INTEGER(np)[0], thinn = INTEGER(thin)[0];
  double prop, acc_prob, post_ratio;
  double *crow, *prow, dpst_lower, dpst_upper;
  SEXP ans, nacc, nex, mc, current;

  nr = 1 + ftrunc(nn/thinn);
  crow = (double *)R_alloc(nnp, sizeof(double));
  prow = (double *)R_alloc(nnp, sizeof(double));
  PROTECT(current = allocVector(REALSXP, nnp));
  PROTECT(nacc = allocVector(REALSXP, nnp));
  PROTECT(nex = allocVector(REALSXP, nnp));
  PROTECT(mc = allocVector(REALSXP, nr * nnp));
  PROTECT(ans = allocVector(VECSXP, 3));

  for(i=0;i<nnp;i++) {
    prow[i] = REAL(init)[i];
    REAL(mc)[i] = REAL(init)[i];
    REAL(nex)[i] = REAL(nacc)[i] = 0.0;
  }

  GetRNGstate();
  for(i=0;i<nn;i++) {
    for(j=0;j<nnp;j++) {     

      prop = rnorm(prow[j], REAL(propsd)[j]);
      
      for(k=0;k<nnp;k++) {
	if (k < j) 
	  REAL(current)[k] = crow[k];

        else
	  REAL(current)[k] = prow[k];

      }      

      defineVar(install("x"), current, env);
      dpst_lower = REAL(eval(f, env))[0];
      
      REAL(current)[j] = prop;
      defineVar(install("x"), current, env);
      dpst_upper = REAL(eval(f, env))[0];
      
      post_ratio = exp(dpst_upper - dpst_lower);
      
      if (dpst_upper == MINF)
        REAL(nex)[j] = REAL(nex)[j] + 1;

      acc_prob = fmin2(1.0, post_ratio);

      if (!R_FINITE(acc_prob)){
	acc_prob = 0.0;
	warning("NaN returned for posterior density");
      }

      if (runif(0, 1) < acc_prob) {
        crow[j] = prop;
	REAL(nacc)[j] = REAL(nacc)[j] + 1;
      }

      else crow[j] = prow[j];
    }

    if( ((i+1) % thinn) == 0){
      l++;
      for(j=0;j<nnp;j++)
	REAL(mc)[l * nnp + j] = crow[j];
    }
    
    for (j=0;j<nnp;j++)
      prow[j] = crow[j];
  }

  PutRNGstate();
  SET_VECTOR_ELT(ans, 0, mc);
  SET_VECTOR_ELT(ans, 1, nacc);
  SET_VECTOR_ELT(ans, 2, nex);
  UNPROTECT(5);
  return(ans);
}
