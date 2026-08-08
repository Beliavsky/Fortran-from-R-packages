#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
  #include <string.h>
  #include <R_ext/Rdynload.h>
  typedef void fun_eval(double x[], int n, double *f, void *data);
  typedef void grad_eval(double x[], int n, double grad[], void *data);
  typedef void callback(double x[], int n, double f, size_t iter, void *data);
  int minimize_nonneg_cg(double x[], int n, double *fun_val,
                         fun_eval *obj_fun, grad_eval *grad_fun, callback *cb, void *data,
                         double tol, size_t maxnfeval, size_t maxiter, size_t *niter, size_t *nfeval,
                         double decr_lnsrch, double lnsrch_const, size_t max_ls,
                         int extra_nonneg_tol, double *buffer_arr, int nthreads, int verbose);
}

typedef struct fholder {
  Rcpp::Function eval_fun;
  Rcpp::Function eval_grad;
  Rcpp::Function docaller;
  Rcpp::List fargs;
} fholder;

void wrapped_feval(double x[], int n, double *f, void *data)
{
  fholder* fdata = (fholder*) data;
  Rcpp::NumericVector temp(fdata->docaller(fdata->eval_fun, fdata->fargs));
  f[0] = temp[0];
}

void wrapped_geval(double x[], int n, double grad[], void *data)
{
  fholder* fdata = (fholder*) data;
  Rcpp::NumericVector temp(fdata->docaller(fdata->eval_grad, fdata->fargs));
  memcpy(grad, &temp[0], sizeof(double) * n);
}

// [[Rcpp::export]]
Rcpp::List run_minimizer(Rcpp::Function feval, Rcpp::Function geval, Rcpp::Function docall, Rcpp::NumericVector x0,
                         Rcpp::List fargs, double tol, size_t maxnfeval, size_t maxiter, double decr_lnsrch,
                         double lnsrch_const, size_t max_ls, int extra_nonneg_tol, int nthreads, int verbose)
{
  int n = x0.size();
  Rcpp::NumericVector buffer_arr(n * 4);
  double *x = x0.begin();
  fargs["x"] = x0;
  // fholder fpass = {
  //   .eval_fun = feval,
  //   .eval_grad = geval,
  //   .docaller = docall,
  //   .fargs = fargs,
  // };
  fholder fpass = {
    feval,
    geval,
    docall,
    fargs,
  };

  size_t slot_niter, slot_nfeval;
  double funval;

  int term_type = minimize_nonneg_cg(x, (int) n, &funval,
                     wrapped_feval, wrapped_geval, NULL, (void*) &fpass,
                     tol, maxnfeval, maxiter, &slot_niter, &slot_nfeval,
                     decr_lnsrch, lnsrch_const, max_ls,
                     extra_nonneg_tol, buffer_arr.begin(), nthreads, verbose);
  Rcpp::List res;
  res["x"] = x0;
  res["nfeval"] = slot_nfeval;
  res["niter"] = slot_niter;
  res["fun"] = funval;
  if (term_type == 0){
    res["term"] = "tol_achieved";
  } else if (term_type == 1){
    res["term"] = "stop_maxnfeval";
  } else if (term_type == 2){
    res["term"] = "stop_maxiter";
  } else {
    res["term"] = "error";
  }
  return res;
}

