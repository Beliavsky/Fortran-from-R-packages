# mlr-fortran v0.1.0

Modern Fortran/FPM translation of the computational framework portions of the
R package **mlr 2.19.3**.

`mlr` is primarily an orchestration package: 150 `RLearner_*` files adapt
algorithms implemented by other R packages. Those external algorithms are not
source code contained in mlr and are therefore not reimplemented here under
mlr's name. The Fortran port focuses on computations that mlr itself performs,
plus a small useful native learner set.

## Implemented computational areas

- Regression measures: SSE, MSE, RMSE, median SE, SAE, MAE, median AE, R2,
  explained variance, RRSE, RAE, MAPE, MSLE/RMSLE, Spearman rho, Kendall tau.
- Classification measures: confusion matrix, MMCE/accuracy, BER/BAC, Cohen
  kappa and weighted kappa, binary confusion counts/rates, precision/NPV/FDR,
  MCC, F1, geometric mean, AUC, Brier/scaled Brier, log loss and multiclass
  Brier score.
- Aggregations: mean, sample SD, RMSE aggregation and `.632` bootstrap.
- Resampling: holdout, k-fold CV, repeated k-fold CV, repeated subsampling and
  bootstrap with out-of-bag test indices.
- Preprocessing: standardization, mean/median/constant imputation, clipping,
  one-hot encoding, class downsampling and oversampling.
- SMOTE: nearest-neighbor construction plus mixed numeric/categorical-code
  synthetic sample generation, translating the package's native `smote.c`
  computation.
- Learners implemented directly: featureless regression/classification,
  ordinary/weighted linear regression, binary logistic regression, Lloyd
  k-means, k-NN regression and classification.
- Survival: Cox proportional-hazards fitting/prediction, Harrell concordance
  and Kaplan-Meier through the supplied `survival-fortran-v0.1.0` dependency.
- Generic callback-based resampling/evaluation for regression and
  classification learners.
- Hyperparameter search: grid search and random search over numeric vectors.
- Feature selection: exhaustive, greedy forward and random subset search.
- Binary probability-threshold tuning.

All project source is free-form Fortran with implicit typing and implicit
external procedures disabled.

## Build

```text
fpm build
fpm test
fpm run --example basic_workflow
```

The package has a path dependency on the included
`vendor/survival-fortran-v0.1.0`, which in turn contains its spline dependency.

## Basic example

```fortran
use mlr_kinds, only : dp, i8
use mlr_rng, only : rng_state, rng_seed
use mlr_resampling, only : make_kfold
use mlr_evaluate, only : resample_regression
use mlr_metrics, only : measure_rmse
```

See `example/basic_workflow.f90` for a complete callback-driven 5-fold
regression example.

## Scope boundary

This is a computational translation, not a clone of R's S3/S4/data-frame
framework. Formula parsing, package discovery, XML, `data.table`, plotting,
parallelMap, R learner registration, and the 150 wrappers around third-party
learners are deliberately not recreated. `API_MAPPING.md` gives the detailed
boundary and future-work candidates.

## Licensing

The original mlr package is BSD-2-Clause. The supplied survival Fortran
translation is GPL-2.0-or-later. The combined linked distribution is therefore
provided under GPL-2.0-or-later while mlr-derived source retains its BSD terms.
See `LICENSE`, `LICENSE.upstream`, `UPSTREAM.md`, and the vendored dependency
notices.
