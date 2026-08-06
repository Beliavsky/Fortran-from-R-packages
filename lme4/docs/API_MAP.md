# API coverage map

## Directly represented computational functionality

| R `lme4` area | Fortran API | Notes |
|---|---|---|
| `lmer` / `mkLmerDevfun` / `optimizeLmer` | `fit_lmm` | Dense marginal covariance; ML and REML |
| penalized least-squares predictor module | `fit_lmm_pls` | Woodbury formulation; avoids an `n by n` covariance matrix |
| `glmer` / PIRLS / Laplace fitting | `fit_glmm` | Five built-in family/link combinations |
| custom `family` objects | `family_spec_t`, `fit_glmm_custom` | Procedure callbacks for link, inverse link, derivative, variance, likelihood, and validity |
| `glmer.nb` | `fit_glmer_nb` | Bounded log-scale profile of negative-binomial size |
| `glmer(..., nAGQ > 1)` scalar | `fit_glmm_aghq` | One grouped scalar random coefficient |
| multidimensional `nAGQ` | `fit_glmm_aghq_multidimensional` | One grouped correlated random-effect vector; tensor quadrature |
| `nlmer` numerical core | `fit_nlmm` | Gaussian response, callback nonlinear mean, one grouped random-effect vector |
| `fixef` | `result%beta` | Fixed-effect estimates |
| `ranef` | `result%u`, `random_effects_for_term` | Conditional modes |
| `VarCorr` / structured covariance concepts | `result%varcorr`, covariance constants | Unstructured, diagonal, compound symmetry, AR(1) |
| `sigma` | result scale fields | Gaussian residual SD or GLMM dispersion/size |
| `fitted`, `residuals` | result arrays | Conditional fitted values and response residuals |
| `predict` | prediction routines | Caller supplies numeric design/covariate arrays |
| `simulate` / `simulate.formula` numerical core | simulation routines | Optional deterministic seed |
| `confint` Wald/profile concepts | `wald_confint_*`, `profile_confint_*_beta` | Fixed effects |
| `bootMer` parametric mode | `parametric_bootstrap_*` | Refits simulated responses |
| percentile bootstrap intervals | `bootstrap_percentile_confint` | Rows are statistics, columns are simulations |
| `influence.merMod` group deletion | `influence_*_groups` | DFBETA and Cook distance by grouping level |
| model comparison | `likelihood_ratio_test` | Chi-square likelihood-ratio test |
| `lmList` numerical fitting | `fit_lm_list`, `predict_lm_list` | Explicit group indices and numeric design |
| `GHrule` | `gh_rule` | Standard-normal nodes, weights, log densities |
| covariance conversion | `sdcor2cov`, `cov2sdcor`, relative-factor routines | Typed matrix API |
| `rePCA` | `re_pca`, `covariance_pca` | Symmetric eigendecomposition |
| singular-fit checking | `is_singular` | Eigenvalue-based test |
| optimizer layer | bundled `minqa` | BOBYQA for multidimensional bounded optimization |

## Built-in GLMM families

| Constant | Mean/link | Variance convention |
|---|---|---|
| `family_binomial` | logit | `mu*(1-mu)` with prior weights |
| `family_poisson` | log | `mu` |
| `family_gamma` | log | `dispersion*mu**2` |
| `family_inverse_gaussian` | log | `dispersion*mu**3` |
| `family_negative_binomial` | log | `mu + mu**2/size` |

Custom family constructors include `gaussian_identity_family`, `binomial_probit_family`, `binomial_cloglog_family`, and `quasipoisson_log_family`. Callers may also populate a `family_spec_t` with their own procedures.

## Changed implementation or narrower scope

- No formula parser is involved. The caller constructs all fixed/random designs and group indices.
- `fit_lmm` uses the dense marginal covariance; `fit_lmm_pls` uses dense Woodbury/PLS matrices. Neither is a sparse CHOLMOD implementation.
- General GLMM fitting uses dense PIRLS and a first-order Laplace objective.
- Multidimensional AGHQ supports one grouping factor and one vector-valued random-effect term. Multiple crossed quadrature blocks are not implemented.
- Tensor AGHQ scales exponentially with random-effect dimension and enforces `max_nodes`.
- `fit_nlmm` currently supports a Gaussian response and one grouping factor, with derivatives of the user mean calculated numerically.
- Fixed-effect profiling drops one design column and refits nuisance parameters. Variance-component profile intervals are not provided.
- Parametric bootstrap is implemented; semiparametric residual and case bootstrap modes are not.
- Prediction requires caller-created numeric arrays and does not build contrasts or handle unseen factor levels automatically.

## R infrastructure intentionally not translated

- Formula parsing, model frames, factors, contrasts, missing-value actions, and `reformulas`
- S3/S4/reference classes and methods
- Sparse Eigen/CHOLMOD storage, symbolic analysis, sparse downdates, and very-large-model workflows
- Data-frame/tibble reconstruction, printing, summaries, plotting, lattice/grid output, and documentation helpers
- `allFit` across external optimizer ecosystems, NLopt wrappers, and optimizer-specific diagnostics
- MCMC sampling and deprecated `mcmcsamp` infrastructure
- Multidimensional AGHQ across several crossed/nested random-effect terms
- Specialized nonlinear response modules beyond Gaussian Laplace fitting
