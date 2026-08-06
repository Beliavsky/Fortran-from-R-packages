# tsmarch-fortran

Modern Fortran translation of the computational core of the R package
`tsmarch` 1.0.3.

The library implements feasible multivariate ARCH workflows for dynamic
conditional correlation (DCC), copula-GARCH, and generalized orthogonal GARCH
(GO-GARCH).  It is designed for FPM and also includes direct GNU Fortran build
scripts for environments where FPM is not installed.

## Implemented numerical scope

- DCC and asymmetric DCC recursion with arbitrary alpha, gamma, and beta orders
- Constant-correlation multivariate Gaussian and Student-t likelihoods
- Separate univariate GARCH marginal estimation through the bundled `tsgarch`
  translation
- Gaussian and Student copula-GARCH probability transforms, filtering,
  estimation, and simulation-based forecasts
- FastICA (symmetric and deflation variants) and a deterministic RADICAL-style
  entropy-rotation implementation
- GO-GARCH estimation and simulation-based forecasting
- Conditional covariance/correlation arrays, coskewness, cokurtosis, and
  portfolio moment aggregation
- EWMA and Ledoit-Wolf covariance estimates, positive-semidefinite repair,
  multivariate normal/Student densities and simulation
- Engle-Sheppard constant-correlation diagnostic
- Simulation and Gaussian value-at-risk/expected-shortfall calculations
- Deterministic density convolution and `dfft`/`pfft`/`qfft` evaluation
- Numerical Hessian, score, covariance, standard-error, AIC, and BIC output for
  DCC fits

## Deliberately omitted

The following are R presentation or object-system features rather than isolated
numerical kernels and are not reproduced:

- S3 dispatch, formulas, model frames, `xts`/`zoo` date indexing, and data-table
  containers
- plotting, `flextable` output, print/summary formatting, and package datasets
- future/parallel orchestration and R-specific sandwich-method dispatch
- direct binary interfaces to Rcpp, RcppArmadillo, RcppParallel, and RcppBessel

Typed Fortran derived types replace R lists and S3 objects.

## Build

With GNU Fortran:

```sh
make check
make optimized
make demo
```

With FPM:

```sh
fpm test
fpm run --example demo_tsmarch
```

The checked build uses strict Fortran 2018 conformance, warnings as errors,
bounds checking and backtraces.  The optimized build uses
`-O3` while retaining strict diagnostics.

## Minimal example

```fortran
use ghyp_kinds, only : dp
use tsgarch, only : garch_spec, fit_options
use tsmarch

type(dcc_spec) :: multivariate_spec
type(garch_spec) :: marginal_spec
type(fit_options) :: options
type(dcc_fit) :: fit

marginal_spec%model = 'garch'
marginal_spec%distribution = 'std'
multivariate_spec%distribution = 'mvt'
multivariate_spec%alpha_order = 1
multivariate_spec%beta_order = 1
options%compute_inference = .true.

fit = estimate_dcc(returns, marginal_spec, multivariate_spec, options)
```

See `example/demo_tsmarch.f90` for a complete deterministic workflow.

## Numerical equivalence

This is a computational translation, not a binary reimplementation of the R
runtime.  The following substitutions are intentional:

- bounded derivative-free optimization and finite-difference inference from the
  GPL-2-compatible `tsdistributions`/`tsgarch` translations replace `nloptr`,
  `Rsolnp`, and `numDeriv` calls;
- FastICA and RADICAL use self-contained deterministic numerical algorithms;
- the portable convolution implementation uses a direct DFT rather than an
  external FFT library;
- time-series indexes and names are not part of numerical outputs.

Consequently, fitted values should be numerically comparable but need not be
bit-for-bit identical to R, especially for stochastic optimization, ICA sign or
ordering conventions, and simulation output.

## Licensing

`tsmarch` declares GPL-2.  This translation is distributed under
GPL-2.0-only so that it remains compatible with the upstream package and the
vendored computational dependencies.

The user-supplied `nloptr` Fortran translation is retained only under
`provenance/dependencies/`.  It is LGPL-3.0-or-later and is not compiled or
linked, because LGPL-3 code cannot be statically combined into a GPL-2-only
work under those exact license choices.

See `DEPENDENCY_NOTES.md`, `PORTING_NOTES.md`, and `licenses/GPL-2.0.txt`.
