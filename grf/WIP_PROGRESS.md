# Parity progress

This checkpoint records the computational parity work completed after the initial GRF translation. It intentionally preserves in-progress work rather than claiming a finished compatibility release.

## Completed before this checkpoint

- Automatic OOB regression-forest nuisance estimates for causal and instrumental forests, including `Z.hat` for IV.
- Automatic OOB multi-response nuisance estimates for `lm_forest` and multi-treatment paths.
- Grouped GRF-style half-sample / little-bag training via `ci_group_size`, with complete-group tree-count rounding.
- Objective-Bayes-debiased grouped prediction variances for regression, causal, IV, causal-survival, and local-linear predictions.
- Causal-survival event/censoring/treatment nuisance forests, IPCW augmentation, doubly robust numerator/denominator moments, scores, and grouped variances.
- Row-specific treatment-variance causal scores and auxiliary binary-instrument compliance forests for IV orthogonal scores.
- Boosted-regression automatic stopping with OOB debiased error and reuse of tuned first-stage controls.
- Regression hyperparameter tuning with deterministic mini-forests and a self-contained Gaussian-process surrogate.
- Local-linear OOB ridge selection and optional residual-splitting training.

## Added in the current checkpoint

- Cluster-first tree sampling with canonical cluster IDs.
- Cluster-level honesty partitioning and cluster-aware OOB status.
- `equalize_cluster_weights` support that samples the same number of observations from each selected cluster and rejects simultaneously supplied explicit sample weights.
- Cluster-robust ATE, calibration, and best-linear-projection standard errors, plus cluster-level RATE half-sample bootstrap.
- Categorical multi-arm detection for non-baseline one-hot treatment contrasts.
- OOB probability-forest estimation of all treatment-arm propensities for categorical multi-arm forests.
- Arm-specific AIPW score construction relative to the baseline arm; the prior general matrix-treatment score remains available for noncategorical multivariate treatments.
- Survival training accepts a custom strictly increasing failure-time grid and stores relabeled time indices.
- Survival splitting uses an event-balance constraint and an exact log-rank statistic; `fast_logrank` selects a streamlined accelerated approximation.
- Survival prediction supports Kaplan-Meier survival, Nelson-Aalen survival/cumulative hazard, arbitrary evaluation grids, and one evaluation time per query row.
- `stabilize_splits` is plumbed through `lm_forest` and `multi_arm_causal_forest` with upstream-compatible defaults.
- Stabilized causal/IV/causal-survival splitting now uses parent-node treatment/instrument means and weighted parent variances, enforces child counts on both sides of the parent mean, applies the alpha-scaled child-variance constraint, and uses the variance-based GRF imbalance penalty.
- `stabilize_splits = .false.` now falls back to the generic regression-style split constraints instead of retaining causal/IV balance checks.
- Stabilized LM/multi-arm splitting enforces the corresponding parent-mean and parent-variance constraints independently for every centered treatment dimension, while retaining the multi-causal count-based imbalance penalty.
- All twelve upstream `generate_causal_data` benchmark selectors are implemented: `simple`, `aw1`, `aw2`, `aw3`, `aw3reverse`, `ai1`, `ai2`, `kunzel`, `nw1`, `nw2`, `nw3`, and `nw4`.
- Causal simulations now apply the upstream `sigma.m`, `sigma.tau`, and `sigma.noise` scaling semantics through optional `sigma_m`, `sigma_tau`, and `sigma_noise` arguments.
- All six upstream `generate_causal_survival_data` selectors are implemented: `simple1` and `type1` through `type5`.
- Causal-survival simulations support selector-specific horizons, caller-supplied predictors, latent-normal AR(1) predictor correlation, deterministic Monte Carlo targets, and treatment-effect signs.

## Still in progress

- Family-specific automatic tuning outside the regression family remains incomplete.
- Some specialized survival/causal-survival inference and edge-case behavior remains streamlined.
- Exact multi-arm covariance/inference variants, alternate HC/TMLE/subset inference options, and the complete upstream RATE/AUTOC/QINI weighting/bootstrap variants remain partial.
- Exact R RNG streams and C++ threading behavior are not translated.

The package therefore remains `substantial`, not `complete`.
