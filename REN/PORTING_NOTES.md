# Porting notes

## Deliberate fixes and clarified behavior

1. The upstream rolling analysis allocates an `LW` method column but never fills
   it. This port fills that slot with the already-exported `po.covShrink`
   computation.
2. Standard maximum drawdown is computed from the running wealth peak. The R
   code compares every observation with the full-sample maximum, which is not
   the conventional definition.
3. `prepare_data` accepts a typed integer `YYYYMMDD` date vector and a separate
   return matrix. It normalizes the first observed month to month 1 and removes
   every return column containing the missing sentinel.
4. The source's one-month weight timing is preserved by default through
   `analysis_options%legacy_weight_timing = .true.`. Set it to `.false.` to use
   weights estimated from months `i-6:i-1` during month `i`.
5. The R parallel loops are serial in this portable implementation. Randomized
   ensembles are reproducible through explicit seeds.

## Numerical replacements

- `quadprog::solve.QP` with only the budget equality uses a pseudoinverse-based
  global minimum-variance solution.
- Long-only minimum variance uses accelerated projected gradient on the unit
  simplex after positive-definite regularization.
- `Matrix::nearPD` is replaced by corpcor's eigenvalue-flooring
  `make_positive_definite`.
- `glmnet::cv.glmnet` for Gaussian models is replaced by a local standardized
  cyclic coordinate-descent path with deterministic K-fold mean-squared-error
  selection. Only the Gaussian subset used by REN is included.
- Canonical correlation for clustering is obtained from principal-angle
  singular values of centered column spaces.
