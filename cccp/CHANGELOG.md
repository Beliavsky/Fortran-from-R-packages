# Changelog

## 0.3.3-fortran.1

- Translated the complete computational surface of cccp 0.3-3.
- Added LP, QP, nonlinear convex, SOCP, SDP, geometric-programming, L1, and
  risk-parity interfaces.
- Replaced Rcpp/RcppArmadillo and the CVXOPT-derived primal-dual implementation
  with a self-contained primal log-barrier/Newton method using BLAS/LAPACK.
- Added phase-I strict-interior recovery for linear, conic, and nonlinear
  constraints.
- Added cone Jordan algebra helpers and approximate primal/dual result fields.
- Added FPM packaging, examples, strict tests, and reproducible build scripts.
