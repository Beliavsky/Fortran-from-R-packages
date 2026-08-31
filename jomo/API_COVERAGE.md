# API coverage

This document maps the numerical surface of upstream **jomo 2.7-6** to the
Fortran implementation. R formula/data-frame/S3/printing/plotting code is out
of scope by design.

| Upstream area | Fortran coverage | Notes |
| --- | --- | --- |
| `jomo1con`, `jomo1cat`, `jomo1mix` | Implemented | Shared by `jomo1_mixed_mcmc`; continuous/categorical convenience wrappers included. |
| `jomo1rancon`, `jomo1rancat`, `jomo1ranmix` | Implemented | `jomo1ran_mixed_mcmc` supports general random-effect design columns and full random-effect covariance. |
| `jomo1ran*hr` | Implemented major kernel | Cluster-specific covariance, hierarchical Wishart scale, and sampled `a` Metropolis update are present. |
| `jomo2com*` | Implemented | Mixed level-1 and level-2 latent-normal data and joint random-effect/level-2 covariance. |
| `jomo2hr*` | Implemented | Direct two-level heterogeneous driver combines joint random/level-2 covariance with cluster-specific level-1 covariance and sampled `a`. |
| `jomo1smc*` proposal/rebuild loop | Implemented | `smc_level1_sweep` uses the upstream `omega(k,k)/10` symmetric random walk for substantive predictors and exact conditional-normal Gibbs draws for auxiliaries. |
| `jomo1ransmc*` / `jomo1ranhrsmc*` compatibility | Implemented | The same SMC driver rebuilds both fixed and random substantive designs; full random-intercept/slope likelihoods and random-effect covariance updates are supported. Cluster-specific level-1 covariances are accepted directly. |
| `jomo2smc*` / `jomo2hrsmc*` compatibility | Implemented | `smc_level2_sweep` updates missing cluster-level predictors against the joint-model density plus the substantive likelihood for all affected rows. |
| SMC powers/interactions/categorical expansion | Implemented | `smc_design_spec` represents the numerical content of upstream `submod`, `ordersub`, and `submodran` without R formula parsing. Continuous powers, categorical `K-1` indicators, and arbitrary interactions are rebuilt after each proposal. |
| Gaussian substantive model (`lm`/`lmer`) | Implemented | Fixed effects, optional random effects, random covariance, Gaussian residual variance, and missing substantive outcomes. |
| Binary probit substantive model (`glm`/`glmer`) | Implemented | Probit latent response, fixed/random effects, random covariance, and missing substantive outcomes. |
| Ordinal probit substantive model (`polr`/`clmm`) | Implemented | Latent response, ordered thresholds, fixed/random effects, random covariance, and missing substantive outcomes. |
| Cox substantive model (`coxph`) | Implemented | Ordered-risk-set partial likelihood and upstream-style coordinate Newton coefficient update. Upstream SMC Cox has no substantive random effects. |
| One complete SMC-compatible iteration | Implemented | `smc_compatible_iteration` performs level-1 proposals, optional level-2 proposals, full design rebuilding, and substantive-parameter updating for a supplied current joint-model mean/covariance state. |
| MCMC-chain variants | Implemented where useful | Single-level optional beta/covariance chains; other joint-model engines return final state/posterior means rather than reproducing R list layouts. |
| Latent categorical augmentation | Implemented | Uses jomo's multivariate latent-normal identification and rejection constraints. |
| Missing-data conditional MVN | Implemented | Explicit logical masks replace R `NA` tests. |
| Wishart/MVN support | Implemented | New free-form Fortran implementation with deterministic RNG state. |
| R formula parsing and `lme4`/`survival`/`ordinal` object extraction | Not translated | Interface-specific; numeric design specifications/matrices are supplied directly. |
| Printing, summaries, plotting, S3 methods, tibbles | Not translated | R-specific or non-computational. |

## SMC parity details

The compatibility step follows the native C algorithms rather than merely
conditioning the joint model after a standard imputation. For a missing latent
coordinate that enters the substantive fixed or random design, the proposal is
symmetric normal with mean equal to the current value and variance
`covariance(k,k)/10`. The Metropolis target combines the corresponding
joint-model Gaussian quadratic form with the substantive likelihood. For a
missing coordinate not used by the substantive model, the driver performs the
exact univariate conditional-normal Gibbs draw. Categorical covariates are
represented by their latent-normal blocks and are decoded/re-expanded whenever
a proposal changes their category.

For Cox models the full ordered-risk-set partial likelihood is recomputed for a
proposal, preserving the fact that changing one row's covariates can change
risk-set denominators for earlier events. The implementation favors a shared,
clear numerical driver over the upstream source's repeated optimized code
blocks; it may therefore do more work per proposal while targeting the same
SMC density.

## Remaining meaningful parity targets

1. Add optional complete-chain storage to all multilevel, heterogeneous, and
   two-level result types, matching the information content of the upstream
   `.MCMCchain` variants.
2. Add higher-level helpers that convert a formula-like description directly
   into `smc_design_spec`. This is convenience/interface parity; polynomial,
   interaction, categorical, level-2, and random-effect numerical expansions
   are already implemented.
3. Add specialized incremental Cox risk-set caches and row-local likelihood
   shortcuts for performance. The current implementation recomputes the exact
   substantive likelihood after a substantive-predictor proposal.

No listed routine silently falls back to an unrelated model.
