#include "header.h"

void gevlik(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns){

  //It computes the log-likelihood for the GEV
  double iscale = 1 / *scale, ishape = 1 / *shape;

  if( (*scale <= 0) || (*shape < -1)) {
    *dns = -1e6;
    return;
  }

  double dummy;

  if (fabs(*shape) <= 1e-16){
    for (int i=0;i<*n;i++){
      if (!ISNA(data[i])){
	dummy = (data[i] - *loc) * iscale;
	*dns += log(iscale) - dummy - exp(-dummy);
      }
    }
  }

  else{
    for(int i=0;i<*n;i++){

      if (!ISNA(data[i])){
	dummy = 1 + *shape * (data[i] - *loc) * iscale;

	if (dummy <= 0) {
	  *dns = -1e6;
	  return;
	}

	*dns += log(iscale) - R_pow(dummy, -ishape) -
	  (ishape + 1) * log(dummy);
      }
    }
  }

  return;
}

void gpdlik(double *exceed, int *n, double *thresh, double *scale,
	    double *shape, double *dns){
  //It computes the log-likelihood for the GPD
  double iscale = 1 / *scale, ishape = 1 / *shape;

  if ((*scale <= 0) || (*shape < -1)) {
    *dns = -1e6;
    return;
  }

  double dummy;

  if (fabs(*shape) <= 1e-16){
    for (int i=0;i<*n;i++){
      dummy = (exceed[i] - *thresh) * iscale;

      if (dummy <= 0){
	*dns = -1e6;
	return;
      }

      *dns += log(iscale) - dummy;
    }
  }

  else{
    for (int i=0;i<*n;i++) {
      dummy = (exceed[i] - *thresh) * iscale;

      if (dummy <= 0) {
	*dns = -1e6;
	return;
      }

      dummy = 1 + *shape * dummy;

      if (dummy <= 0) {
	*dns = -1e6;
	return;
      }

      *dns += log(iscale) - (ishape + 1) * log(dummy);
    }
  }

  return;
}
