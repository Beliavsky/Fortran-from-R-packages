# Translation coverage

## Included

- exponential negative log-likelihood, gradient, prediction, penalty, and
  proximal solver;
- Gamma fixed-shape negative log-likelihood and gradient;
- joint Gamma coefficient/shape maximum likelihood;
- FISTA or nonaccelerated proximal-gradient iteration with backtracking;
- L1 and elastic-net proximal operators and objective values;
- lambda maximum and logarithmic grids;
- fixed-lambda fits and warm-start regularization paths;
- repeated deterministic balanced K-fold cross-validation;
- mean squared-error helpers from `FISTA.R`;
- compatibility entry points corresponding to the exported Rcpp functions.

## Replaced or corrected

- the zero-vector `FitGlmFixed` C++ stub is replaced by working fixed-lambda
  and path APIs;
- the nonimplemented C++ Gamma likelihood stubs are replaced by the working R
  Gamma formulation;
- invalid upstream CV index arithmetic is replaced by balanced folds;
- positive Gamma shape is enforced through a log parameter;
- the correct elastic-net proximal map is the default, with the C++ expression
  retained as an option.

## Omitted

- Rcpp/RcppEigen registration and conversion code;
- R documentation machinery and package-loading code;
- interface-demo arithmetic in `MyClass`;
- examples requiring external R data packages or `RPEIF` objects.

The original package tree is preserved under `upstream/` for auditability.
