# ranger — modern Fortran computational translation

This directory is a modern free-form Fortran translation of the computational algorithms in the R package **ranger 0.18.0**, a fast random-forest implementation by Marvin N. Wright and contributors.

The public module is `ranger`. The translation is designed to live as a top-level sibling of `rfortran-core` in the `Fortran-from-R-packages` repository. It imports the shared real kind `dp` from `r_kinds`; the package does not define a second floating-point kind and does not vendor Rcpp, RcppEigen, Matrix, BLAS/LAPACK, or another translated R package.

## Build

From the repository layout

```text
Fortran-from-R-packages/
  rfortran-core/
  ranger/
```

run:

```text
cd ranger
fpm build
fpm test
fpm run --example classification_example
fpm run --example survival_example
```

## Main API

The public module exports typed forest and option objects together with:

- `fit_ranger_classification`, `predict_ranger_classification`
- `fit_ranger_probability`, `predict_ranger_probability`
- `fit_ranger_regression`, `predict_ranger_regression`
- `predict_ranger_quantiles`
- `fit_ranger_survival`, `predict_ranger_survival`
- Janitza and Altmann variable-importance p-value helpers
- infinitesimal-jackknife prediction variance
- hierarchical shrinkage for regression and probability forests
- case-specific terminal-node weights, tree sizes, and variable-use counts

The `ranger_options` type controls tree count, `mtry`, node and bucket sizes, depth, replacement and sample fraction, split rule, ExtraTrees candidate count, split-selection weights, always-split variables, regularization, importance mode, in-bag storage, quantile-regression node-value preparation, holdout evaluation, `oob_error`, and missing-value routing.

The four forest families retain ranger-specific computational distinctions rather than being aliases for one generic tree:

- classification: Gini, Hellinger, and ExtraTrees splits;
- probability: class-probability leaves with classification-style splitting;
- regression: variance, ExtraTrees, Maxstat, Beta, and Poisson splits;
- survival: log-rank, ExtraTrees, AUC/C-index-style, and Maxstat splits with Nelson–Aalen terminal cumulative hazards. `fit_ranger_survival` accepts either an explicit `time_interest` vector (sorted and deduplicated) or a scalar count of approximately equally spaced observed event times.

Numeric and integer-coded categorical predictors are supported. Unordered factors can be treated as ordered or searched by partitions. Missing numeric/factor predictor values can use learned left/right routing.

Maxstat splitting includes ranger's native rank/logrank scoring, Lau92/Lau94 p-value corrections, Benjamini-Hochberg adjustment, and the survival single-cut unadjusted-normal special case. See `API_COVERAGE.md` for the remaining implementation-level adaptations.

OOB error computation can be disabled with `options%oob_error=.false.` without disabling permutation importance. When OOB computation is enabled but an observation receives no OOB trees, regression and probability OOB predictions follow ranger and are NaN; survival CHF rows remain zero. Classification uses integer OOB labels, so label 0 plus `oob_count==0` represents ranger's no-prediction case.

## Design notes

The R formula/data-frame/S3/Rcpp presentation layer is intentionally omitted. Native callers supply dense numeric predictor arrays and integer class/factor codes directly. See `API_COVERAGE.md` for exact translated, adapted, and omitted areas.

For maintainability, ranger's C++ low-level cache/memory-mode and small-Q/large-Q split-search micro-optimizations are not copied literally. Candidate splits are evaluated directly in modern Fortran using the same statistical criteria. This favors clarity and reproducibility over reproducing C++-specific storage optimizations.

## Licensing and provenance

The upstream R package declares GPL-3. Its C++ core explicitly states that the C++ core is MIT licensed and carries copyright (c) 2014-2018 Marvin N. Wright. This translation is distributed as part of the GPL-3 package while retaining the MIT core notice and attribution in `LICENSE.MIT-CORE` and `NOTICE.md`.

See `NOTICE.md`, `UPSTREAM.md`, `upstream/`, and `API_COVERAGE.md` for detailed provenance.
