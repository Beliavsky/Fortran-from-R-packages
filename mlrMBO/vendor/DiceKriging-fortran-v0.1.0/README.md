# DiceKriging-fortran

A modern Fortran/FPM translation of the computational core of the R package
**DiceKriging 1.6.1** (2025-10-20).

The goal is numerical functionality, not emulation of R's S4/S3, formula,
data-frame, graphics, or printing infrastructure. The implementation is
self-contained and has no R, BLAS/LAPACK, `rgenoud`, `foreach`, or other runtime
dependency.

## Implemented computational functionality

- `km`-style Gaussian-process/kriging model fitting with explicit trend matrices.
- Tensor-product Gaussian, exponential, Matern 3/2, Matern 5/2, and
  power-exponential covariance models.
- Isotropic covariance models.
- Nonstationary input scaling with piecewise-linear inverse-range fields.
- Process variance, nugget, and heteroskedastic observation-noise handling.
- Concentrated/profile maximum likelihood (MLE).
- Penalized MLE with the SCAD penalty.
- Leave-one-out (LOO) fitting, LOO criterion, analytic LOO gradient, and LOO
  predictions.
- Analytic covariance and likelihood derivatives for stationary kernels;
  finite-difference covariance derivatives for the scaling parameterization.
- Simple and universal kriging prediction, conditional covariance, bias
  correction, and 95% prediction intervals.
- Conditional and unconditional Gaussian-process simulation.
- Sequential model update with fixed covariance parameters or full covariance
  re-estimation.
- Response-only updates for already-existing design points.
- Cross-validation by arbitrary fold labels.
- Covariance-vector spatial derivatives and trend derivatives.
- Branin, six-hump camelback, Goldstein-Price, Hartman-3, and Hartman-6 test
  functions.

In particular, the Fortran equivalent of the behavior needed by
`KrigInv::CovReEstimate=TRUE` is:

```fortran
call km_update(model, newx, newy, newf, cov_reestimate=.true.)
```

This appends the observations and re-runs covariance-parameter estimation,
rather than keeping the original ranges fixed.

## Basic use

```fortran
use dicekriging

type(km_model) :: model
type(km_prediction) :: pred
real(dp), allocatable :: f(:,:)

call trend_constant(x, f)
call km_fit(model, x, y, f, 'matern5_2')
call km_predict(model, xnew, fnew, 'UK', pred, se_compute=.true.)
```

`x` is `n x d`; `f` is the corresponding `n x p` trend/design matrix. Explicit
matrices replace R formulas and `model.matrix`. Helpers are supplied for the
most common trends: constant, linear, linear with two-factor interactions, and
quadratic.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example fit_branin
```

The release was also built directly with gfortran 14.2.0 using Fortran 2018,
runtime checking, warnings-as-errors, and implicit-interface checking. See
`VALIDATION.md`.

## Important compatibility notes

- DiceKriging uses R's `optim`/L-BFGS-B and optionally `rgenoud`. This port uses
  a native bounded multistart BFGS optimizer. The fitted numerical optima are
  validated against upstream DiceKriging regression values, but optimizer
  trajectories are not expected to be identical.
- R formula parsing, S4/S3 classes/method dispatch, data-frame/name checking,
  printing, summaries, and plotting are intentionally omitted.
- `covUser` accepts an arbitrary R function as a covariance kernel. Arbitrary
  R callbacks are language-binding infrastructure and are not represented by
  this self-contained library. The built-in DiceKriging covariance families
  are implemented natively.
- `matern5_2add0` is referenced by `covStruct.create.R`, but its
  `covAdditive0` class is not supplied by DiceKriging 1.6.1 itself; it is not
  implemented here.
- Parallel `foreach`/`doParallel` execution is not reproduced. Independent
  starts/folds can be parallelized by a calling Fortran application if desired.

## License

The supplied DiceKriging package declares `GPL-2 | GPL-3`. This translation is
provided under the same choice. The complete GPL version 2 and GPL version 3
texts are included as `LICENSE-GPL-2` and `LICENSE-GPL-3`. Upstream package
metadata is retained under `upstream/`.

This is an independent source translation and is not an official release of
the DiceKriging authors.
