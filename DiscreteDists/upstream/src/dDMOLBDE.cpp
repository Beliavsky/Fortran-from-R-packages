#include <Rcpp.h>
#include "shared.h"
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;

inline double logpmf_dDMOLBE(double x, double mu, double sigma,
                            bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) )
    return x+mu+sigma;
#endif
  if (mu < 0.0 || sigma < 0.0 || x < 0.0 || !isInteger(x, false)) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(x) || x < 0.0)
    return R_NegInf;

  double k1 = (1+x/mu) * std::exp(-x/mu);
  double k2 = (1+(x+1)/mu) * std::exp(-(x+1)/mu);
  double p1 = std::log(sigma) + std::log(k1 - k2);
  double p2 = -std::log(1-(1-sigma)*k1) - std::log(1-(1-sigma)*k2);
  double res = p1 + p2;

  return res;
}

// [[Rcpp::export]]
NumericVector cpp_dDMOLBE(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& log_prob = false
) {

  if (std::min({x.length(), mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);

  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_dDMOLBE(GETV(x, i),
                          GETV(mu, i),
                          GETV(sigma, i),
                          throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);

  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}
