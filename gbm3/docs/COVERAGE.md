# Computational coverage

## Translated

### Boosting/tree engine

- Initial function estimates.
- Working responses/pseudo-residuals.
- Stochastic bagging by observation ID.
- Random predictor subsampling.
- Best-first tree growth.
- Continuous/categorical splits.
- Missing-value child nodes.
- Monotonicity constraints.
- Distribution-specific terminal-node updates.
- Shrinkage updates.
- Training and validation deviance histories.
- Out-of-bag improvement histories.
- Stored raw fitted values.
- Continuation training (`gbm_more` computational workflow).
- Full prediction and response-scale transformations.
- Per-tree contributions and staged multi-tree-count prediction.
- Partial dependence (`gbm_plot` computational kernel).
- Friedman split-based relative influence.
- Permutation relative influence.
- Friedman's H interaction statistic.
- Tree/node inspection helpers.

### Cross-validation and performance

- Automatic k-fold CV orchestration.
- User-supplied fold IDs.
- Repeated observation IDs kept in one fold.
- Pairwise ranking groups kept in one fold.
- Optional Bernoulli class stratification.
- Fold-size-weighted CV error curves.
- Out-of-fold fitted values at the selected iteration.
- Training and validation/test best-iteration selectors.
- Raw cumulative-OOB best-iteration selector.

### Distribution families

- Gaussian.
- Bernoulli.
- Poisson.
- Gamma.
- Laplace/absolute loss.
- Student-t loss and location-M terminal fitting.
- Quantile loss.
- AdaBoost exponential loss.
- Huberized hinge loss.
- Tweedie loss.
- Cox proportional hazards, including right-censored and counting-process data.
- Pairwise/LambdaMART ranking.

### Pairwise ranking measures

- NDCG.
- Concordance.
- Mean average precision (MAP).
- Mean reciprocal rank (MRR).
- Rank cutoffs (`pairwise_max_rank`).

### Survival

- Efron and Breslow tie handling in boosting likelihood calculations.
- Cox baseline-hazard increments and cumulative baseline hazard.
- Cox continuation, CV, and permutation-importance extension.

## Intentionally omitted R/interface functionality

- Formula parsing and model-frame construction.
- Data-frame/factor conversion and R attributes/classes.
- S3 printing/summary methods.
- Plotting, calibration plots, rug plots, and performance plots.
- Rcpp/SEXP entry points and registration.
- R parallel-backend orchestration and progress output.

## Remaining native API gaps

- The exact LOESS-smoothed OOB selector used by `gbm.perf(method="OOB")` is
  not reproduced because it depends on R's LOESS implementation. The Fortran
  API exposes the unsmoothed calculation explicitly as `oob_raw`.
- Portable, versioned save/load serialization has not yet been defined. The
  public derived types can be inspected directly, but compiler-specific memory
  dumps are intentionally not presented as a file format.

## Exact-parity qualifications

This is a source translation, not an R runtime emulator. With the same data
and deterministic tree decisions, formulas target upstream behavior. Exact
bit-for-bit identity is not claimed where:

1. R's RNG affects bagging, feature order, CV fold randomization, permutation
   order, or rank tie-breaking.
2. Cox risk-set sums occur in a different floating-point order.
3. OpenMP execution order differs from the serial Fortran implementation.
4. R-specific preprocessing would encode or reorder data before entering the
   C++ core.
5. The explicitly named `oob_raw` selector is used instead of R LOESS.
