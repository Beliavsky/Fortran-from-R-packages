# Changelog

## 1.1.4-fortran.1

- translated all computational smoothing, bandwidth-selection, derivative, confidence-bound, ARMA, forecasting, bootstrap, and rolling-backtest algorithms
- recovered and encoded the hidden `sysdata.rda` algorithm lookup tables
- replaced Rcpp/Armadillo smoothing and forecasting kernels with modern Fortran
- added a self-contained conditional Gaussian ARMA estimator
- added typed result structures and status codes
- added FPM application, examples, four test programs, and independent fixed references
- retained original source and GPL-3.0-only licensing for provenance
