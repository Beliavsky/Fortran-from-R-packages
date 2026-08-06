#include <Rcpp.h>
using namespace Rcpp;

// Disable floating-point contraction (fused multiply-add) so the recursion is
// evaluated with the same separate multiply/add roundings as the plain R loop,
// making the C++ output bit-for-bit identical to it.
#pragma STDC FP_CONTRACT OFF

// GARCH(1,1) residual recursion for span_simulate().
// Given unit-variance innovations z, returns eps_t = sigma_t * z_t with
// sigma^2_{t+1} = omega + alpha * eps_t^2 + beta * sigma^2_t, initialised at the
// unconditional variance omega / (1 - alpha - beta). The innovations are drawn
// in R, so the RNG stream is unchanged and the output is identical to the plain
// R recursion; only the sequential loop moves to C++.
//
// [[Rcpp::export]]
NumericVector garch_filter(NumericVector z, double omega, double alpha, double beta) {
  const R_xlen_t n = z.size();
  NumericVector eps(n);
  double s2 = omega / (1.0 - alpha - beta);
  for (R_xlen_t t = 0; t < n; ++t) {
    double e = std::sqrt(s2) * z[t];
    eps[t] = e;
    s2 = omega + alpha * e * e + beta * s2;
  }
  return eps;
}
