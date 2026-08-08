#include "evaluate.h"
#include <limits>
using namespace arma;

thread_local Rcpp::EvalBase *fev = NULL;     // pointer to abstract base class
thread_local Rcpp::EvalBase *gev = NULL;     // pointer to abstract base class

#include "n1qn1.h"

extern "C" void n1qn1_ (S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
                        int mode[], int niter[], int nsim[], int imp[], double zm[], int izs[], float rzs[], double dzs[], int []);

thread_local unsigned int nq1n1c_calls = 0, nq1n1c_grads = 0;
thread_local int nq1n1c_fprint = 0;
static void fwrap(int *ind, int *n, double *x, double *f, double *g, int *ti, float *tr, double *td, int *id)
{
  int i;
  Rcpp::NumericVector par(*n), ret(*n);
  std::copy(&x[0], &x[0]+*n, &par[0]);

  if (*ind==2 || *ind==4) {
    nq1n1c_calls++;
    ret = fev->eval(par);
    if (nq1n1c_fprint){
      Rprintf("%3d:%#14.8g:", nq1n1c_calls, ret[0]);
      for (i = 0; i < *n; i++) Rprintf(" %#8g", x[i]);
      Rprintf("\n");
    }
    *f = ret[0];
  }
  if (*ind==3 || *ind==4) {
    nq1n1c_grads++;
    ret = gev->eval(par);
    std::copy(&ret[0],&ret[0]+*n,&g[0]);
    // for (i = 0; i < *n; i++) g[i] = ret[i];
  }
}

uvec lowerTri(mat H, bool diag = false){
  unsigned int d = H.n_rows;
  mat o(d, d, fill::ones);
  if (!diag){
    return find(trimatl(o,-1));
  } else {
    return find(trimatl(o));
  }
}


RcppExport SEXP n1qn1_wrap(
           SEXP fSEXP, SEXP gSEXP, SEXP rhoSEXP, SEXP xSEXP, SEXP epsSEXP,
           SEXP nSEXP, SEXP modeSEXP, SEXP niterSEXP, SEXP nsimSEXP, SEXP impSEXP,
           SEXP nzmSEXP, SEXP zmSEXP, SEXP fprint_sexp) {
  BEGIN_RCPP
    nq1n1c_calls=0;
  nq1n1c_grads=0;
  nq1n1c_fprint = INTEGER(fprint_sexp)[0];

  // Validate inputs before allocating fev/gev: Rf_error longjmps out,
  // so any allocations made before the error check would be leaked.
  int n, mode, niter, nsim, imp, nzm;
  n = INTEGER(nSEXP)[0];
  mode = INTEGER(modeSEXP)[0];
  niter = INTEGER(niterSEXP)[0];
  nsim = INTEGER(nsimSEXP)[0];
  imp = INTEGER(impSEXP)[0];
  nzm = INTEGER(nzmSEXP)[0];

  if (n <= 0) Rf_error("n must be positive");
  if (nzm <= 0) Rf_error("nzm must be positive (check for integer overflow in nzm calculation)");
  // Prevent signed int overflow in triangular index arithmetic n*(n+1)/2
  // used in n1qn1_all.c:144 and below.  sqrt(INT_MAX) ~ 46340.95.
  if ((double)n * (n + 1) > (double)INT_MAX) {
    Rf_error("n is too large: n*(n+1) overflows 32-bit integer (maximum n is ~46340)");
  }

  delete fev;
  if (TYPEOF(fSEXP) == EXTPTRSXP){
    fev = new Rcpp::EvalCompiled(fSEXP, rhoSEXP); // xptr
  } else {
    fev = new Rcpp::EvalStandard(fSEXP, rhoSEXP); // Standard evaulation
  }
  delete gev;
  if (TYPEOF(gSEXP) == EXTPTRSXP){
    gev = new Rcpp::EvalCompiled(gSEXP, rhoSEXP); // xptr
  } else {
    gev = new Rcpp::EvalStandard(gSEXP, rhoSEXP); // Standard evaulation
  }

  double f, eps;
  double *x = new double[n];
  double *g = new double[n];
  double *var = new double[n];
  double *zm = new double[nzm];

  int izs[1]; float rzs[1]; double dzs[1];
  std::copy(&(REAL(xSEXP)[0]),&(REAL(xSEXP)[0])+n, &x[0]);
  std::copy(&(REAL(zmSEXP)[0]),&(REAL(zmSEXP)[0])+nzm, &zm[0]);
  eps = REAL(epsSEXP)[0];
  std::fill(&var[0], &var[0]+n, 0.1);
  int id = 0;

  n1qn1_(fwrap,&n,x,&f,g,var,&eps,
         &mode,&niter,&nsim,&imp,zm,izs,rzs,dzs, &id);

  Rcpp::NumericVector par(n);
  std::copy(&x[0],&x[0]+n,&par[0]);
  // for (i=0; i<n; i++) par[i] = x[i];
  Rcpp::NumericVector hess(nzm);
  // On input this is hessian
  // On output this is H = LDL'
  // Triangular matrix is paramterized by column instead of row.
  mat L = eye(n,n);
  mat D = mat(n,n,fill::zeros);
  mat H = mat(n,n);
  arma::uword zmLen = (arma::uword)n * (n + 1) / 2;
  vec zmV(zmLen);
  std::copy(&zm[0], &zm[0] + zmLen, zmV.begin());
  H.elem(lowerTri(H,true)) = zmV;
  if (n == 1) H = D;
  else{
    L.elem(lowerTri(H,false)) = H.elem(lowerTri(H,0));
    D.diag() = H.diag();
    H = L*D*L.t();
  }
  // Hessian -> c.hess
  vec hessV = H.elem(lowerTri(H,true));
  std::copy(hessV.begin(),hessV.end(),hess.begin());

  delete[] x;
  delete[] g;
  delete[] var;
  delete[] zm;

  return Rcpp::List::create(Rcpp::Named("value") = f,
                            Rcpp::Named("par") = par,
                            // Rcpp::Named("L") = L,
                            // Rcpp::Named("D") = D,
                            Rcpp::Named("H") = H,
                            // Rcpp::Named("zm")=zms,
                            Rcpp::Named("c.hess") = hess,
			    Rcpp::Named("n.fn") = nq1n1c_calls,
			    Rcpp::Named("n.gr") = nq1n1c_grads);
  END_RCPP
 }
