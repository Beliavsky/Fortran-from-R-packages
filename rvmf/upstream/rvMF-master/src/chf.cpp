#include <Rcpp.h>
#include "mpreal.h"

// This function is from the scModels package.
// kummer fn: Taylor series method
double kummer_taylor(double x, double a, double b) {
  mpfr::mpreal aj = 1;
  mpfr::mpreal sj = aj;
  mpfr::mpreal tol = 1e-6, err = 1.0, j = 0.0;
  mpfr::mpreal aj1 = 0.0, sj1 = 0.0;
  mpfr::mpreal x_mp = mpfr::mpreal(x), a_mp = mpfr::mpreal(a), b_mp = mpfr::mpreal(b);
  while (err > tol) {
    aj1 = aj*(a+j)*x/((b+j)*(j+1));
    sj1 = sj + aj1;
    aj = aj1;
    sj = sj1;
    err = mpfr::abs(aj1);
    j = j+1;
  }
  return mpfr::log(sj).toDouble();
}

// computes log(M(a = d1/2, b = d1, z = 2*k))
// [[Rcpp::export(.chf)]]
double chf(double k, double d1) {
  return kummer_taylor(2.0*k, d1/2.0, d1);
}
