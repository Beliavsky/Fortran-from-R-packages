#define USE_FC_LEN_T

#include "R.h"
#include "Rinternals.h"
#include "Rmath.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

typedef double (scaleFnPtr)(int, double*, double*, double*, double*);
typedef double (rcovFnPtr)(int, double*, double*, scaleFnPtr*, double*, double*, double*);


double my_mad(int n, double *x, double *dwork1, double *dwork2, double *mu);
double my_median(int n, double *x);
double gk(int n, double *x, double *y, scaleFnPtr *scalefn, double *dwork1, double *dwork2, double *dwork3);
double qc(int n, double *x, double *y, scaleFnPtr *scalefn, double *dwork1, double *dwork2, double *dwork3);
double dsum(int n, double* x, int incx, double* wrkn);
double scaleTau2(int n, double *x, double *dwork1, double *dwork2, double *mu);


SEXP covOPW(SEXP SX, SEXP Siter, SEXP SscaleFun, SEXP SrcovFun)
{
  char CHARA = 'A', CHARL = 'L', CHARN = 'N', CHART = 'T', CHARV = 'V';
  double *X = NULL, *Z = NULL, *ZCOPY = NULL, *U = NULL, **A = NULL, *d = NULL;
  double *dwork1 = NULL, *dwork2 = NULL, *dwork3 = NULL, *diagT = NULL, *offdiagT = NULL;
  double *tau = NULL, *gamma = NULL, *cov = NULL, *covcopy = NULL, *center = NULL, *dist = NULL;
  double mu = 0.0, alpha = 0.0, DZERO = 0.0, DONE = 1.0;
  int n = 0, p = 0, np = 0, pp = 0, iter = -1, i = 0, j = 0, k = 0, info = 0, lwork = 0;
  int liwork = 0, IONE = 1;
  int *isuppz = NULL, *iwork = NULL;
	SEXP Sans = R_NilValue, Scov = R_NilValue, Scenter = R_NilValue;
	SEXP Sdist = R_NilValue, Sdim = R_NilValue, Snames = R_NilValue;
  scaleFnPtr *scalefn = NULL;
  rcovFnPtr *rcovfn = NULL;

  if(strncmp(CHAR(asChar(SscaleFun)), "s_mad", 5) == 0)
    scalefn = &my_mad;
  else if(strncmp(CHAR(asChar(SscaleFun)), "scaleTau2", 9) == 0)
    scalefn = &scaleTau2;
  else
    error("unable to set scale function pointer in C function covOPW");

  if(strncmp(CHAR(asChar(SrcovFun)), "gk", 2) == 0) 
    rcovfn = &gk;
  else if(strncmp(CHAR(asChar(SrcovFun)), "qc", 2) == 0) 
    rcovfn = &qc;
  else
    error("unable to set rcov function pointer in C function covOPW");

  if(!isMatrix(SX))
    error("first argument to C function covOPW is not a matrix");

  PROTECT(Sdim = getAttrib(SX, R_DimSymbol));
  n = INTEGER(Sdim)[0];
  p = INTEGER(Sdim)[1];
  Sdim = R_NilValue;
  UNPROTECT(1);

  np = n*p;
  pp = p*p;
  lwork = 18*p;
  liwork = 10*p;

  iter = INTEGER(Siter)[0];

  X = REAL(SX);
  Z = (double*) R_alloc((size_t) np, sizeof(double));
  F77_CALL(dcopy)(&np, X, &IONE, Z, &IONE);

  ZCOPY = (double*) R_alloc((size_t) np, sizeof(double));
  U = (double*) R_alloc((size_t) (p*(p+1))/2, sizeof(double));
  covcopy = (double*) R_alloc((size_t) pp, sizeof(double));

  A = (double**) R_alloc((size_t) iter, sizeof(double*));
  for(k = 0; k < iter; k++)
    A[k] = (double*) R_alloc((size_t) pp, sizeof(double));

  d = (double*) R_alloc((size_t) p, sizeof(double));
  dwork1 = (double*) R_alloc((size_t) n, sizeof(double));
  dwork2 = (double*) R_alloc((size_t) n, sizeof(double));
  dwork3 = (double*) R_alloc((size_t) imax2(n, lwork), sizeof(double));
  diagT = (double*) R_alloc((size_t) p, sizeof(double));
  offdiagT = (double*) R_alloc((size_t) p, sizeof(double));
  tau = (double*) R_alloc((size_t) (p-1), sizeof(double));
  gamma = (double*) R_alloc((size_t) p, sizeof(double));
  isuppz = (int*) R_alloc((size_t) (2*p), sizeof(int));
  iwork = (int*) R_alloc((size_t) liwork, sizeof(int));



  for(k = 0; k < iter; k++) {

    for(j = 0; j < p; j++) {
      d[j] = scalefn(n, Z+j*n, dwork1, dwork2, &mu);
      //this can be handled better
      if(fabs(d[j]) < 1e-12)
        error("column with zero scale encountered in C function covOPW");
      alpha = 1.0 / d[j];
      F77_CALL(dscal)(&n, &alpha, Z+j*n, &IONE);
    }

    for(i = 0; i < p; i++)
      U[i+((2*p-i-1)*i)/2] = 1.0;

    for(i = 1; i < p; i++)
      for(j = 0; j < i; j++)
        U[i+((2*p-j-1)*j)/2] = rcovfn(n, Z+i*n, Z+j*n, scalefn, dwork1, dwork2, dwork3);

/*	Rprintf("\n %d (%d): %d, %d, %d, %d %d %d \n", k, iter, n, p, np, pp, lwork, liwork);
*/
    F77_CALL(dsptrd)(&CHARL, &p, U, diagT, offdiagT, tau, &info FCONE);
    F77_CALL(dstegr)(&CHARV, &CHARA, &p, diagT, offdiagT, &mu, &mu, &i,
                     &i, &DZERO, &j, gamma, A[k], &p, isuppz, dwork3,
                     &lwork, iwork, &liwork, &info FCONE FCONE);
    F77_CALL(dopmtr)(&CHARL, &CHARL, &CHARN, &p, &p, U, tau, A[k], &p,
                     dwork2, &info FCONE FCONE FCONE);

    for(j = 0; j < p/2; j++)
      F77_CALL(dswap)(&p, A[k]+j*p, &IONE, A[k]+p*(p-j-1), &IONE);

    F77_CALL(dcopy)(&np, Z, &IONE, ZCOPY, &IONE);
    F77_CALL(dgemm)(&CHARN, &CHARN, &n, &p, &p, &DONE, ZCOPY, &n, A[k],
                    &p, &DZERO, Z, &n FCONE FCONE);

    for(i = 0; i < p; i++)
      for(j = 0; j < p; j++)
        A[k][i+j*p] = d[i] * A[k][i+j*p];
  }

  PROTECT(Scov = allocMatrix(REALSXP, p, p));
  PROTECT(Scenter = allocVector(REALSXP, p));
  PROTECT(Sdist = allocVector(REALSXP, n));
  cov = REAL(Scov);
  center = REAL(Scenter);
  dist = REAL(Sdist);

  for(j = 0; j < p; j++) {
    gamma[j] = scalefn(n, Z+j*n, dwork1, dwork2, &mu);
    for(i = 0; i < p; i++)
      cov[i+j*p] = i == j ? gamma[j] * gamma[j] : 0.0;
    center[j] = mu;
  }

  for(i = 0; i < n; i++) {
    for(j = 0; j < p; j++)
      Z[i+j*n] = R_pow_di(((Z[i+j*n] - center[j]) / gamma[j]), 2);
    dist[i] = F77_CALL(dasum)(&p, Z+i, &n);
  }

  for(k = iter-1; k >= 0; k--) {
    F77_CALL(dcopy)(&pp, cov, &IONE, covcopy, &IONE);
    F77_CALL(dgemm)(&CHARN, &CHARN, &p, &p, &p, &DONE, A[k], &p, covcopy,
                    &p, &DZERO, cov, &p FCONE FCONE);
    F77_CALL(dcopy)(&pp, cov, &IONE, covcopy, &IONE);
    F77_CALL(dgemm)(&CHARN, &CHART, &p, &p, &p, &DONE, covcopy, &p, A[k],
                    &p, &DZERO, cov, &p FCONE FCONE);
    F77_CALL(dcopy)(&p, center, &IONE, gamma, &IONE);
    F77_CALL(dgemv)(&CHARN, &p, &p, &DONE, A[k], &p, gamma, &IONE,
                    &DZERO, center, &IONE FCONE);
  }

  PROTECT(Sans = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(Sans, 0, Scenter);
  SET_VECTOR_ELT(Sans, 1, Scov);
  SET_VECTOR_ELT(Sans, 2, Sdist);
	PROTECT(Snames = allocVector(STRSXP, 3));
	SET_STRING_ELT(Snames, 0, mkChar("center"));
	SET_STRING_ELT(Snames, 1, mkChar("cov"));
	SET_STRING_ELT(Snames, 2, mkChar("distances"));
	setAttrib(Sans, R_NamesSymbol, Snames);
  UNPROTECT(5);

  return Sans;
}


double my_median(int n, double *x)
{
  int k = -1;

  if(n%2) {
    k = (n-1)/2;
    rPsort(x, n, k);
    return(x[k]);
  }

  k = n/2;
  rPsort(x, n, k);
  rPsort(x, k, k-1);
  return((x[k-1] + x[k]) / 2.0);
}


double my_mad(int n, double *x, double *dwork1, double *dwork2, double *mu)
{
  const int IONE = 1;
  int i = 0;
  F77_CALL(dcopy)(&n, x, &IONE, dwork1, &IONE);
  *mu = my_median(n, dwork1);
  for(i = 0; i < n; i++)
    dwork1[i] = fabs(dwork1[i] - *mu);
  return my_median(n, dwork1) * 1.4826;
}


double gk(int n, double *x, double *y, scaleFnPtr *scalefn, double *dwork1, double *dwork2, double *dwork3)
{
  const int IONE = 1;
  double plus = 0.0, minus = 0.0, mu = 0.0;
  const double DONE = 1.0, DNEGONE = -1.0;

  F77_CALL(dcopy)(&n, x, &IONE, dwork1, &IONE);
  F77_CALL(daxpy)(&n, &DONE, y, &IONE, dwork1, &IONE);
  plus = scalefn(n, dwork1, dwork2, dwork3, &mu);

  F77_CALL(dcopy)(&n, x, &IONE, dwork1, &IONE);
  F77_CALL(daxpy)(&n, &DNEGONE, y, &IONE, dwork1, &IONE);
  minus = scalefn(n, dwork1, dwork2, dwork3, &mu);

  return (R_pow_di(plus, 2) - R_pow_di(minus, 2)) / 4.0;
}


double scaleTau2(int n, double *x, double *dwork1, double *dwork2, double *mu)
{
  const double C1 = 4.5, C2squared = 9.0;
//  const double C2 = 3.0;
  const double Es2c = 0.9247153921761315;
	double medx = 0.0, sigma0 = 0.0, tmpsum = 0.0;
  int i = 0, IONE = 1;

  F77_CALL(dcopy)(&n, x, &IONE, dwork1, &IONE);
  medx = my_median(n, dwork1);
  for(i = 0; i < n; i++)
    dwork1[i] = fabs(dwork1[i] - medx);
  sigma0 = my_median(n, dwork1);

  F77_CALL(dcopy)(&n, x, &IONE, dwork1, &IONE);
	for(i = 0; i < n; i++) {
		dwork1[i] = fabs(dwork1[i] - medx);
		dwork1[i] = dwork1[i] / (C1 * sigma0);
		dwork2[i] = 1.0 - R_pow_di(dwork1[i], 2);
		dwork2[i] = R_pow_di(((fabs(dwork2[i]) + dwork2[i])/2.0), 2);
  }

  tmpsum = dsum(n, dwork2, 1, dwork1);

  for(i = 0; i < n; i++)
		dwork1[i] = x[i] * dwork2[i];

	*mu = dsum(n, dwork1, 1, dwork2) / tmpsum;

  F77_CALL(dcopy)(&n, x, &IONE, dwork1, &IONE);
  for(i = 0; i < n; i++) {
    dwork2[i] = R_pow_di((dwork1[i] - *mu) / sigma0, 2);
    dwork2[i] = dwork2[i] > C2squared ? C2squared : dwork2[i];
  }

  return sigma0 * sqrt(dsum(n, dwork2, 1, dwork1) / (n*Es2c));
}


double qc(int n, double *x, double *y, scaleFnPtr *scalefn, double *dwork1, double *dwork2, double *dwork3)
{
  double medx = 0.0, medy = 0.0, r = 0.0;
  int IONE = 1, i = 0, onethree = 0, twofour = 0;

  F77_CALL(dcopy)(&n, x, &IONE, dwork1, &IONE);
  medx = my_median(n, dwork1);
  F77_CALL(dcopy)(&n, y, &IONE, dwork1, &IONE);
  medy = my_median(n, dwork1);

  for(i = 0; i < n; i++) {
    if((x[i] > medx && y[i] > medy) || (x[i] < medx && y[i] < medy))
      onethree++;
    else if((x[i] > medx && y[i] < medy) || (x[i] < medx && y[i] > medy))
      twofour++;
  }

  r = ((double) (onethree - twofour)) / ((double) (onethree + twofour));
  return sin(M_PI_2*r);
}


double dsum(int n, double* x, int incx, double* wrkn)
{
  int i = 0;

  if(n == 1)
    return x[0];

  while(i < n/2) {
    wrkn[i] = x[2*incx*i] + x[(2*i+1)*incx];
    i++;
  }

  if(2*i < n)
    wrkn[i-1] = wrkn[i-1] + x[2*incx*i];

  return dsum(i, wrkn, 1, wrkn+i);
}


