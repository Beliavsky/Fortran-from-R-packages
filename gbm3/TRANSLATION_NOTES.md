# Translation notes

## Source target

- Upstream package: `gbm3`
- Upstream version: 3.0.3
- Upstream DESCRIPTION date: 2026-07-15
- Upstream license: GPL (>= 2)
- Translation license: GPL-2.0-or-later

The supplied upstream C++ and R sources are retained in `upstream/`.

## Main source mapping

| Modern Fortran | Principal upstream sources |
| --- | --- |
| `gbm3_tree.f90` | `tree.*`, `node*`, `varsplitter*`, continuous/categorical splitter strategies |
| `gbm3_distributions.f90` | Gaussian, Bernoulli, Poisson, Gamma, Laplace, tdist, quantile, AdaBoost, Huberized, Tweedie sources |
| `gbm3_pairwise.f90` | `pairwise.*` LambdaMART and IR measures |
| `gbm3_cox.f90` | `coxph.*`, censored/counting Cox state implementations |
| `gbm3_core.f90` | `gbm_engine.*`, `gbm_fit.*`, `gbm_more`, prediction and `gbm_plot` computational behavior |
| `gbm3_cv.f90` | `gbm-cross-val.r`, `create-cv-groups.r`, `gbm-cv-errors.r`, `gbm-cv-predict.r` |
| `gbm3_diagnostics.f90` | `gbm-perf.r`, `permutation-relative-influence.r`, `gbm-interactions.r`, `pretty-gbm-tree.r` |
| `gbm3_math.f90` | location-M, sorting, stable elementary numerical helpers |

## Deliberate Fortran adaptations

- C++ virtual distribution classes are represented by integer distribution
  identifiers and `select case` dispatch.
- R/C++ callbacks, Rcpp objects, SEXP values, formula handling, and data-frame
  conversion are absent from the native API.
- Trees use allocatable derived types and 1-based Fortran indices internally.
- Categorical data retains the upstream data coding convention: category
  values are zero-based integers stored in the real predictor matrix.
- The R RNG is replaced by the standard Fortran RNG. Algorithms and sampling
  rules are retained, but equal integer seeds do not imply identical random
  streams between R and Fortran.
- The Cox likelihood is evaluated directly from risk sets rather than using
  the upstream pointer-sorted state machine. It implements the same target
  partial likelihood and Efron/Breslow tie formulas, but floating-point
  summation order can differ.
- The implementation is serial. Upstream OpenMP scheduling/parallel plumbing
  is not required for numerical functionality and is not translated.

## Expanded API behavior

- `gbm_continue` appends trees and preserves old training/validation/OOB
  histories. The model now stores raw fitted values so continuation carries
  the same old fit forward instead of reconstructing it from changed input
  predictors.
- `gbm_predict_staged` evaluates multiple cumulative tree counts in one pass.
  `gbm_predict_trees` returns individual tree contributions and corresponds to
  the computational purpose of upstream `single.tree=TRUE` prediction.
- `gbm_cross_validate` generates folds before fitting the full model, then fits
  each fold, matching the high-level order of upstream CV orchestration.
  Repeated IDs and pairwise groups are kept intact. Bernoulli stratification
  and explicit fold IDs are supported.
- CV errors use the same fold-size weighting pattern as `gbm_cv_errors`.
- `gbm_permutation_importance` uses one shuffled row permutation for each
  variable in turn, restoring the original predictor between variables.
- `gbm_interaction_strength` implements Friedman's H construction from
  centered partial-dependence values over observed unique combinations.
- `gbm_best_iteration(...,"oob_raw")` intentionally does not claim the LOESS
  smoothing performed by R's `gbm.perf(method="OOB")`.

## Tree details retained

- Best-first growth: one best terminal-node split is committed per depth step.
- Each committed split creates left, right, and missing children.
- Split improvement uses upstream weighted residual variance reduction.
- Missing observations do not participate in the minimum left/right child
  count test.
- Missing terminal nodes smaller than `min_num_obs_in_node` inherit the parent
  prediction after the upstream left/right weighted adjustment.
- Categorical levels are ordered by weighted residual mean before prefix split
  search.
- Feature tie behavior favors the first feature in the shuffled candidate
  order.
- Friedman relative influence follows the upstream left/right traversal and
  intentionally does not recurse through the missing branch.

## Numerical validation

The maintained Fortran sources are tested with strict GNU Fortran flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fimplicit-none -fcheck=all
```

Tests cover all simple distribution dispatch paths, Gaussian and Bernoulli
end-to-end fits, categorical/missing trees, partial dependence, LambdaMART,
Cox PH, continuation parity against one-shot fits, staged/tree-wise prediction,
repeated-ID and stratified CV, pairwise/Cox CV and continuation, permutation
importance, best-iteration selection, model inspection, and Friedman's H.
See `test/test_gbm3.f90`.
