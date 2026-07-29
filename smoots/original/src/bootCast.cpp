#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// written by Dominik Schulz; version: 14/04/2020

// Arguments:

// X: The observed time series.
// innov: The residual series that is usually obtained from the 'arima()'
//        function.
// ar: The AR-coefficients ordered from i = 1 to p.
// ma: The MA-coefficients ordered from j = 1 to q.
// mu: The estimated mean of the series 'X'.
// Fi: The demeaned residual series; this series works as the empirical
//     distribution, from which innovations are generated for the simulation.
// h: The forecasting horizon; this arguments represents the maximum future
//    time period n + h to calculate point forecasts and forecasting intervals
//    for. The calculations are done for all future periods n + 1 to n + h.

// Function for point forecasts with horizon h
// [[Rcpp::export]]
arma::vec fcastCpp(const arma::vec& X, const arma::vec& innov,
  const arma::rowvec& ar, const arma::rowvec& ma, const double mu,
  const int h) {
  const int p = ar.size();
  const int q = ma.size();
  const int l = std::max(p, q);
  const int n = X.size();
  arma::vec Xsub(l + h);
  arma::vec isub(l + h - 1);
  isub.zeros();
  Xsub.subvec(0, l - 1) = X.subvec(n - l, n - 1) - mu;
  isub.subvec(0, l - 1) = innov.subvec(n - l, n - 1);
  const arma::rowvec maa = arma::reverse(ma);
  const arma::rowvec arr = arma::reverse(ar);
  for (int i = 0 + l; i < l + h; ++i) {
    Xsub.subvec(i, i) = arr * Xsub.subvec(i - p, i - 1) +
      maa * isub.subvec(i - q, i - 1);
  }
  return Xsub.subvec(l, l + h - 1) + mu;
}

// Function for simulated 'true' forecasts
// [[Rcpp::export]]
arma::vec tfcastCpp(const arma::vec& X, const arma::vec& innov,
  const arma::vec& epsBoot, const arma::rowvec& ar, const arma::rowvec& ma,
  const double mu, const int h) {
  const int p = ar.size();
  const int q = ma.size();
  const int l = std::max(p, q);
  const int n = X.size();
  arma::vec Xsub(l + h);
  arma::vec isub(l + h );
  Xsub.subvec(0, l - 1) = X.subvec(n - l, n - 1) - mu;
  isub.subvec(0, l - 1) = innov.subvec(n - l, n - 1);
  isub.subvec(l, l + h - 1) = epsBoot;
  const arma::rowvec maa = arma::reverse(ma);
  const arma::rowvec arr = arma::reverse(ar);
  for (int i = 0 + l; i < l + h; ++i) {
    Xsub.subvec(i, i) = arr * Xsub.subvec(i - p, i - 1) +
      maa * isub.subvec(i - q, i - 1) + isub(i);
  }
  return Xsub.subvec(l, l + h - 1) + mu;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
