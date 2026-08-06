# RPEGLMEN modern Fortran

This project translates the computational core of the R package `RPEGLMEN`
1.1.4 to modern Fortran and provides an FPM package.

The library fits positive responses with:

- exponential generalized linear models with an elastic-net penalty;
- Gamma generalized linear models with an elastic-net penalty;
- joint unpenalized Gamma coefficient and shape estimation;
- logarithmic lambda grids, regularization paths, and repeated balanced
  K-fold cross-validation.

Rcpp/RcppEigen bindings, R argument dispatch, package examples that require
`RPEIF`, and other R runtime infrastructure are omitted. The package contains
no plotting engine.

## Build

With FPM:

```text
fpm test
fpm run --example periodogram_glm_demo
```

Without FPM:

```text
make MODE=checked test
make MODE=optimized test
make MODE=checked example
```

The checked build uses Fortran 2018, warnings, bounds/runtime checks, and
backtraces. Shell and Windows helper scripts are in `scripts/`.

## Basic use

```fortran
use rpeglmen, only : dp, enet_options, fit_result, glmnet_exp

real(dp) :: a(200, 4), b(200)
type(enet_options) :: options
type(fit_result) :: fit

options%alpha = 0.5_dp
options%num_lambda = 50
options%k_fold = 5
options%k_fold_iter = 3
call glmnet_exp(a, b, fit, options)
```

The selected coefficients are in `fit%coefficients`; the selected penalty is
`fit%selected_lambda`; and `fit%lambda_grid`, `fit%cv_mean`, and `fit%cv_sd`
contain the cross-validation path.

For Gamma responses:

```fortran
use rpeglmen, only : fit_glm_gamma_net
call fit_glm_gamma_net(a, b, fit, options)
```

The unpenalized Gamma MLE is available as `fit_glm_gamma_mle`, while
`glmnet_exp_fixed` and `glm_gamma_net_fixed` fit a specified penalty.

## Important numerical choices

The default implementation uses the mathematically correct elastic-net
proximal map, FISTA acceleration, deterministic balanced folds, and mean
negative log-likelihood for validation. Options are available to reproduce
several upstream choices:

- `source_proximal=.true.` uses the threshold formula from the C++ code;
- `use_fista=.false.` reproduces the C++ solver's nonaccelerated update;
- `cv_metric='source'` uses response-scale prediction error;
- `normalize_gradient=.true.` reproduces the optional normalized gradient.

The upstream C++ fold-size expression performs integer division before
multiplication, producing a zero-length training set for ordinary K-fold
settings. That path is not reproduced; this port always creates valid,
balanced folds. See `PORTING_NOTES.md`.

## Dependencies

The Fortran implementation is self-contained and does not require BLAS,
LAPACK, Eigen, Rcpp, RPEIF, or an external optimization library. The RPEIF
import in the original package is used by examples rather than by the fitted
GLM engines.

## License

The upstream package is GPL (>= 2). This translation is distributed under
GPL-2.0-or-later. Full license texts and the complete upstream source snapshot
are included.
