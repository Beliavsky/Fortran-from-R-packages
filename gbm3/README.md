# gbm3-fortran

A modern free-form Fortran translation of the computational core of the R
package **gbm3 3.0.3**, packaged for the Fortran Package Manager (FPM).

The implementation is GPL-2.0-or-later, matching the upstream `GPL (>= 2)`
license. See `NOTICE.md`, `LICENSE`, and `upstream/`.

## Implemented computational features

- Stochastic gradient boosting with best-first regression trees.
- Continuous and categorical predictors, explicit missing-value branches, and
  monotonicity constraints.
- Observation-ID bagging and random feature subsampling.
- Gaussian, Bernoulli, Poisson, Gamma, Laplace, Student-t, quantile,
  AdaBoost, Huberized hinge, and Tweedie losses.
- Cox proportional-hazards boosting for right-censored `(time,status)` and
  counting-process `(start,stop,status)` responses, with Breslow or Efron ties.
- LambdaMART/pairwise ranking with NDCG, concordance, MAP, and MRR metrics.
- Link-scale and response-scale prediction.
- Per-tree prediction contributions and staged prediction at multiple tree
  counts in one call.
- Continuation training corresponding to upstream `gbm_more`.
- Native k-fold cross-validation, including repeated-ID grouping, pairwise
  group preservation, optional Bernoulli stratification, CV error curves, and
  out-of-fold fitted values.
- Training/test best-iteration selection and an explicitly named raw OOB
  selector (`oob_raw`).
- Friedman split-based relative influence and permutation importance.
- Partial dependence corresponding to the computational behavior of upstream
  `gbm_plot` (without plotting).
- Friedman's H interaction statistic built from partial dependence.
- Tree/node inspection helpers.
- Cox baseline-hazard estimation.
- Training, validation, and out-of-bag improvement histories.
- Stored raw fitted values on the link scale, excluding offsets, matching the
  quantity carried forward by upstream continuation fitting.

## Real kind

All maintained Fortran source uses one named real kind:

```fortran
use gbm3_kinds, only : dp
```

`dp` is defined once from `iso_fortran_env::real64`. Real variables use
`real(dp)` and real constants use the `_dp` suffix.

## Basic example

```fortran
program example
   use gbm3
   implicit none
   integer, parameter :: n = 100
   real(dp) :: x(n, 1), y(n)
   real(dp), allocatable :: pred(:)
   type(gbm_options) :: options
   type(gbm_model) :: model
   integer :: i

   do i = 1, n
      x(i, 1) = -2.0_dp + 4.0_dp * real(i - 1, dp) / real(n - 1, dp)
      y(i) = 1.0_dp + 2.0_dp * x(i, 1)
   end do

   options = gbm_options(distribution=GBM_GAUSSIAN, num_trees=100, &
                         interaction_depth=2, min_num_obs_in_node=5, &
                         shrinkage=0.05_dp, bag_fraction=0.5_dp)
   call gbm_set_seed(1234)
   call gbm_fit(x, y, model, options)
   pred = gbm_predict(model, x)
end program example
```

## Expanded native API

Continuation training appends trees while retaining the existing fitted values
and performance histories:

```fortran
call gbm_continue(model, x, y, additional_trees=50)
```

For Cox PH, pass the same rank-2 survival response accepted by `gbm_fit`. For
pairwise ranking, also pass `group=`.

Staged prediction mirrors the vector `n.trees` behavior of upstream prediction:

```fortran
integer :: counts(3)
real(dp), allocatable :: staged(:, :), contribution(:, :)

counts = [25, 50, 100]
staged = gbm_predict_staged(model, x, counts)
contribution = gbm_predict_trees(model, x)
```

`gbm_predict_trees` returns one column per individual tree contribution.

Native cross-validation returns a `gbm_cv_result` containing the complete CV
error curve, selected iteration, fold IDs, and out-of-fold fitted values:

```fortran
type(gbm_cv_result) :: cv
type(gbm_model) :: full_model

call gbm_cross_validate(x, y, cv, n_folds=5, options=options, &
                        full_model=full_model)
print *, cv%best_iteration
```

If `id=` is supplied, all rows sharing an ID stay in the same fold. Pairwise
ranking folds are assigned by `group=`. Bernoulli fits can request
`stratify=.true.`. A user-supplied `fold_id=` vector is also accepted.

Other diagnostics include:

```fortran
integer :: best
real(dp) :: h
real(dp), allocatable :: importance(:)

best = gbm_best_iteration(model, "train")
call gbm_permutation_importance(model, x, y, importance)
h = gbm_interaction_strength(model, x, [1, 2])
```

`gbm_best_iteration` accepts `train`, `validation`/`test`, or `oob_raw`.
Upstream `gbm.perf(method="OOB")` first smooths the OOB curve with R's LOESS;
that R-specific smoothing step is not claimed by `oob_raw`.

## Response types

For ordinary distributions, call `gbm_fit(x, y, ...)` with a vector response.
For pairwise ranking, also pass an integer `group=` vector; each group's rows
must be contiguous, as expected by the upstream LambdaMART implementation.

For Cox PH, pass a rank-2 response to the same generic `gbm_fit`:

```fortran
real(dp) :: surv(n, 2)       ! columns: time, status
! or
real(dp) :: surv_cp(n, 3)    ! columns: start, stop, status

options%distribution = GBM_COXPH
call gbm_fit(x, surv, model, options)
```

## Categorical predictors

`var_classes(j)=0` means predictor `j` is continuous. A positive value gives
the number of categories. Categorical values in `x` are represented by the
same zero-based integer codes used by the upstream C++ core (`0,1,...,K-1`),
stored as `real(dp)`. IEEE NaN denotes a missing predictor value.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example gaussian_example
fpm run --example advanced_api_example
```

The source is standard free-form Fortran 2018 and has no required BLAS/LAPACK
or R dependency.

## Scope

The maintained Fortran code targets numerical/computational behavior. R
formula/data-frame handling, S3 classes, plotting, console progress reporting,
and R/CPP interface code are intentionally not translated. Portable,
versioned model serialization is not yet part of the native API. See
`docs/COVERAGE.md` for computational coverage and known parity limits.
