# rpart-fortran

A modern Fortran translation of the computational core of **rpart 4.1-27**,
the R package for recursive partitioning and regression trees.

This project preserves the upstream GPL licensing and keeps an unmodified copy
of the source package under `upstream/` for provenance and parity review.

## Implemented functionality

- ANOVA/regression trees.
- Classification trees with Gini or information impurity.
- User priors and asymmetric classification loss matrices.
- Poisson trees with empirical-Bayes shrinkage and deviance/square-root error.
- Exponential-survival trees for right-censored and start/stop data.
- Continuous, unordered categorical, and ordered-as-numeric predictors.
- Missing predictors, competitor splits, surrogate splits, and all three
  `usesurrogate` policies.
- Variable-specific split costs.
- Cost-complexity values, CP tables, pruning, and cross-validation errors.
- Cross-validated predictions corresponding to `xpred.rpart`, including the
  full terminal-node response when requested.
- Variable importance from primary and surrogate splits.
- Prediction of values, classes, class probabilities, terminal node IDs,
  full node responses, and node paths.
- `rpart.exp` time-axis transformation for survival data.

The actual R formula/data-frame/S3 layer and plotting/interactive utilities are
not part of the native Fortran API.

## Precision policy

All maintained Fortran uses one public real kind:

```fortran
use rpart, only : dp
```

`dp` is defined exactly once as `real64` in `rpart_kinds`. Maintained source
uses `real(dp)` and `_dp` real constants throughout.

## Build with FPM

```text
fpm build
fpm test
fpm run --example basic_example
fpm run --example survival_example
```

The project has no external library dependency.

## Basic use

```fortran
program demo
   use rpart
   implicit none

   real(dp) :: x(6,1), y(6), pred(6)
   type(rpart_control) :: control
   type(rpart_model) :: model
   integer :: stat

   x(:,1) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
   y = [1.0_dp,1.0_dp,1.0_dp,5.0_dp,5.0_dp,5.0_dp]

   control = rpart_make_control(minsplit=2, minbucket=1, cp=0.0_dp, &
                                xval=0, stat=stat)
   if (stat /= 0) error stop "invalid control"

   call rpart_fit_regression(x, y, model, control=control, stat=stat)
   if (stat /= 0) error stop "fit failed"
   call rpart_predict_values(model, x, pred)
end program demo
```

Use `rpart_make_control` when changing `minsplit` or `minbucket`: it reproduces
R's coupled defaults (`minbucket=round(minsplit/3)` and, when only
`minbucket` is supplied, `minsplit=3*minbucket`). Direct construction of an
`rpart_control` object cannot automatically infer that one field was changed.

## Main API

Fitting:

- `rpart_fit_regression`
- `rpart_fit_classification`
- `rpart_fit_poisson`
- `rpart_fit_survival`
- `rpart_fit_survival_startstop`

Prediction:

- `rpart_predict_values`
- `rpart_predict_class`
- `rpart_predict_proba`
- `rpart_predict_full`
- `rpart_predict_where`
- `rpart_predict_one`
- `rpart_node_path`

Cross-validation prediction:

- `rpart_default_xpred_cp`
- `rpart_xpred_regression`
- `rpart_xpred_classification`
- `rpart_xpred_poisson`
- `rpart_xpred_survival`
- `rpart_xpred_survival_startstop`
- `rpart_xpred_full`

Tree utilities:

- `prune_model`
- `compute_variable_importance`
- `count_nodes`
- `count_splits`
- `rpart_exp_transform_right`
- `rpart_exp_transform_startstop`

Categorical predictors are supplied through the optional integer `ncat(:)`
argument. `ncat(j)=0` means continuous/ordered numeric; a positive value gives
the number of unordered categories, encoded as integer-valued reals 1..ncat.
Non-finite predictor values are treated as missing.

## Cross-validation reproducibility

The R package uses R's RNG when it creates fold labels. This port has its own
deterministic RNG, so a seed does not reproduce R's random grouping bit for
bit. For exact fold control, pass `xgroups=` to fitting or explicit `groups`
to the `rpart_xpred_*` routines. The fold-tree CP rescaling and method-specific
loss calculation follow the upstream algorithm.

## Validation

The test suite contains deterministic checks for:

- regression, classification, and Poisson splits and predictions;
- categorical predictors and surrogate routing;
- upstream rpart's hard-core priors/loss fixture, including its published
  split-improvement values;
- Gini and information classification splits;
- right-censored and counting-process survival transformations/fits;
- cost-sensitive split selection;
- CP pruning and zero-risk trees;
- CP-table cross-validation and `xpred` fold predictions.

See `TRANSLATION_NOTES.md` for the mapping from upstream C/R components and
known interface differences.

## License

Upstream `rpart` declares `GPL-2 | GPL-3`. This translation is distributed
under the same choice. See `LICENSE-GPL-2`, `LICENSE-GPL-3`, and `NOTICE.md`.
