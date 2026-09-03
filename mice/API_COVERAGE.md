# API coverage

This file describes **numerical parity**, not R-language interface parity.
The Fortran API works on explicit arrays rather than R formulas, factors,
S3 objects, data frames, quosures, or `mids` objects.

## Implemented computational kernels

### Core FCS and continuous imputers

- Multi-chain, multi-iteration fully conditional specification driver with a
  target-by-predictor matrix, deterministic seed, observed-cell preservation,
  random observed-donor initialization, and chain means/variances.
- `mean`, `sample`.
- Bayesian normal (`norm`), `norm.nob`, `norm.predict`, and bootstrap normal.
- Predictive mean matching with match types 0, 1, and 2.
- Current `matchindex` and legacy `matcher` donor-selection algorithms.
- `midastouch`, including bootstrap-frequency regression, automatic kappa,
  leave-one-out donor models, distance-powered probabilities, and effective
  donor-pool-size output.
- `mpmm`: multivariate PMM using the leading canonical response direction and
  whole-row donor copying.
- Random-indicator and quadratic PMM numerical kernels.
- Ridge fallback for rank-deficient normal regression.

### Categorical imputers

- Binary Bayesian logistic regression with White-Daniel-Royston augmentation.
- Bootstrap logistic regression.
- Nominal multinomial logistic regression with augmentation and posterior
  probability category draws.
- Dedicated proportional-odds cumulative-logit fitting and prediction for
  `polr`, including ordered threshold parameterization and augmented fitting.
- LDA imputation with observed category priors, category means, pooled
  within-category covariance, posterior probabilities, and category draws.
  Singular covariance handling uses a small documented ridge rather than
  reproducing MASS's R-level rank-dropping machinery.

### MNAR / NARFCS numerical kernels

- `mnar.norm` additive sensitivity offsets combined with the same Bayesian
  normal identifiable model as `norm`.
- `mnar.logreg` additive sensitivity offsets combined with the same augmented
  Bayesian logistic identifiable model as `logreg`.
- The Fortran API accepts the unidentifiable design matrix and delta vector
  directly. R's `ums` string parser and `umx` name resolution are intentionally
  omitted as language/interface machinery.

### Two-level variables

- `2lonly.mean`, including repair of partially missing cluster-constant numeric
  variables and NaN for clusters with no observed response.
- `2lonly.norm` and `2lonly.pmm`: cluster consistency checks, cluster-level
  aggregation, imputation, and expansion back to row level.
- Full heterogeneous `2l.norm` Gibbs sampler from Kasim-Raudenbush as used by
  upstream mice: random intercept option, arbitrary random-slope columns,
  cluster-specific residual precisions, population random-effect mean,
  inverse-Wishart random-effect precision, hierarchical `sigma2.0`/`theta`
  updates, and posterior-predictive imputations.

### Pooling and tests

- `pool.scalar` Rubin 1987 and Reiter 2003 rules.
- Barnard-Rubin finite complete-data degrees of freedom.
- Multivariate Rubin mean/within/between/total covariance pooling.
- Pooled Wald quadratic form, the reusable numerical core behind D1-style
  model tests when caller-supplied coefficient contrasts are available.
- Meng-Rubin D3 statistic from caller-supplied fitted and restricted deviances.

### Missing-data diagnostics

- Explicit IEEE-NaN missingness masks using `ieee_is_nan`.
- `md.pairs` four pairwise count matrices.
- `md.pattern` unique patterns, frequencies, margins, and column ordering.
- `flux`: pobs, influx, outflux, average inbound/outbound, and FICO.
- `quickpred` using pairwise Pearson correlations and usable-case thresholds.
  R's optional Kendall/Spearman modes and name-based include/exclude arguments
  are interface conveniences and are not included in the numeric API.

### Amputation and native helpers

- MCAR multivariate amputation.
- Continuous RIGHT/LEFT/MID/TAIL missingness curves with the upstream
  10,000-standard-normal shift-calibration construction.
- Discrete odds/score-group amputation.
- Normalized Legendre basis from the native C++ kernel.

## Meaningful remaining computational gaps

These are intentionally documented rather than represented as exact parity:

1. `2l.bin` and `2l.lmer`. Upstream delegates the model fit to `lme4`; this
   translation does not duplicate a mixed-model engine until a compatible
   callable API from the repository's `lme4` translation is verified.
2. `2l.pan` and `jomoImpute` orchestration. These are primarily adapters to
   external imputation engines and R model/data structures; no dependency
   source is vendored here.
3. CART and random-forest/literanger/ranger imputers. No compatible callable
   repository API was established during this parity pass, so approximate
   replacement trees were not substituted for the upstream engines.
4. Lasso imputers. A top-level `glmnet` translation exists in the target
   repository, but this package does not invent an interface before a stable
   callable Fortran API is verified.
5. Complete D1/D2 model-object workflows and D3 model refitting. The reusable
   array-level pooling, Wald, and D3 arithmetic is present; formula/model
   extraction and nested-model refitting remain outside the numeric API.
6. Passive formula evaluation, block/multivariate formulas, custom R methods,
   `post` expressions, `ignore` expression handling, and R call dispatch.

The remaining R plotting, printing, `broom`, tidyverse, S3, `mids`/`mira`
object manipulation, file/export, and parallel orchestration code is outside
the requested computational translation scope.
