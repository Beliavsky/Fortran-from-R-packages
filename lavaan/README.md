# lavaan-fortran

Modern Fortran/FPM numerical core inspired by `lavaan` 0.7-2.

Version 0.7.0 extends the v0.6 engine with mixed continuous/ordinal conditional-x statistics, missing-response random-coefficient likelihoods, automatic/name-aware MIIV marker rewriting, SAM joint block covariance and second-order bias propagation, and explicit Satorra-Bentler nested-model difference tests.

The library is matrix-native. It does not attempt to emulate lavaan syntax parsing, S4 classes, model frames, printing, or plotting.

## Core SEM representation

- `ram_model` with RAM A, S and optional intercept/mean vector
- `ram_free_map` for mapping free parameters into RAM matrices
- `ram_sigma`, `ram_mu`
- `ram_from_lisrel`
- `standardized_ram`
- CFA, path, latent-regression, growth, and general recursive/non-recursive SEMs whenever `I-A` is nonsingular

## Estimation

- covariance/mean ML
- GLS
- ULS
- WLS
- DWLS
- raw-data ML
- FIML with arbitrary IEEE-NaN missingness patterns
- multigroup covariance ML and multigroup FIML
- nonlinear equality/inequality constraints and bounds
- BFGS optimization
- numerical Hessians/SEs using vendored `numDeriv-fortran`

## EFA and rotation

`efa_fit_cov` supports:

- principal-axis factor extraction (`PAF`)
- concentrated Gaussian ML factor extraction (`ML`)
- orthogonal or oblique rotation through the vendored `GPArotation-fortran` engine
- varimax, quartimax, geomin, oblimin, Crawford-Ferguson, target, simplimax, bifactor and other criteria supported by that dependency
- uniquenesses, communalities, factor correlations, reproduced covariance, and residual covariance

## SAM-style two-stage estimation

- `sam_fit_cov`
- `sam_fit_data`
- `sam_fix_measurement`
- `sam_propagate_uncertainty`
- `sam_yuan_chan_test` for the single-group/global Yuan-Chan residualized-Gamma scaling calculation
- `sam_continuous_gamma` for continuous ADF moment covariance
- `sam_browne_unbiased_gamma` for Browne's finite-sample unbiased covariance-moment Gamma

The first-stage measurement parameters are estimated and then inserted as fixed quantities into a second-stage structural RAM model. v0.4 additionally finite-differences the stage-2 estimates with respect to stage-1 parameters and adds the delta-method term `J V_stage1 J'` to the conditional stage-2 covariance. Both conditional and propagated covariance/SE results are retained. v0.5 also exposes the matrix-level Yuan-Chan global SAM correction, using `D = I - Delta_gamma P`, the residualized `Gamma_tilde = D Gamma D'`, and the structural-only Satorra-Bentler projection.

## MIIV / 2SLS

- `miiv_2sls` for raw-data equations
- `miiv_2sls_cov` for covariance/moment equations
- first-stage F diagnostics
- Sargan over-identification statistic
- `ram_miiv_candidates` for model-implied instrument screening using the RAM disturbance covariance condition
- `ram_miiv_equations` for automatic equation discovery from nonzero RAM paths, descendant exclusion, relevance screening, and equation-level identification flags
- `ram_miiv_marker_equations` for scaling-indicator/marker rewriting of latent structural equations into observed proxy equations, including transformed coefficients and instrument rescreening against the composite proxy disturbance
- `miiv_estimate_uls`, `miiv_estimate_gls`, `miiv_estimate_2rls`, `miiv_estimate_rls`
- matching `miiv_jacobian_*` routines for the covariance-statistic linearization used by lavaan MIIV variance calculations
- `miiv_vcov_from_gamma` for `H Gamma H'/N` covariance propagation

## Categorical estimation

- thresholds and polychoric correlations
- `categorical_wls_statistics`
  - one statistic vector containing all thresholds followed by polychoric correlations
  - delete-one jackknife asymptotic Gamma matrix
- `categorical_wls_statistics_analytic`
  - the same threshold/correlation statistic vector
  - estimating-function/influence-matrix Gamma without leave-one-out refits
  - threshold uncertainty and its effect on polychoric scores are propagated jointly
- `muthen1984_ordinal` / `categorical_wls_statistics_muthen`
  - explicit stage-1 ordinal threshold scores and stage-2 polychoric scores
  - Muthen block-triangular `A11/A21/A22` bread
  - analytic bivariate-normal boundary and correlation derivatives
  - `N * acov(s)` Gamma, full WLS weights, and DWLS diagonal weights
- full WLS and diagonal DWLS weights
- `fit_ram_pml_ordinal`
- `fit_ram_pml_mixed`
  - continuous-continuous Gaussian pairs
  - ordinal-ordinal rectangle probabilities
  - continuous-ordinal conditional-normal pair likelihoods

The bivariate normal probabilities use vendored `pbivnorm-fortran`.

## Marginal maximum likelihood

`fit_mml_ordinal_factor` implements normal-ogive ordinal factor MML using tensor Gauss-Hermite quadrature for one to three latent factors. Marginal thresholds are fixed at empirical probit estimates, factor variances are standardized, and free loading patterns are supplied explicitly.

v0.4 adds `fit_mml_mixed_factor` / `mml_mixed_loglik` for mixed continuous/ordinal outcomes, arbitrary fixed latent means/covariances, IEEE-NaN missing values, free loading patterns, and continuous-item intercept/residual-scale estimation. `mml_ordinal_loglik` and `gauss_hermite_normal` remain public.

v0.5 adds `mml_mixed_loglik_qmc` and `fit_mml_mixed_factor_qmc`. They use deterministic Halton-normal integration and support one to eight latent dimensions, providing a practical high-dimensional alternative when tensor Gauss-Hermite quadrature is no longer feasible.

v0.6 adds `mml_mixed_loglik_adaptive` and `fit_mml_mixed_factor_adaptive` for one to three latent dimensions. Each observation gets a posterior mode and Hessian, and Gauss-Hermite nodes are recentered/rescaled to the local Gaussian posterior approximation. The resulting importance-quadrature form is exact for Gaussian conditional models and is much more efficient than nonadaptive tensor quadrature when the posterior is concentrated.

## Two-level SEM

Two paths are available:

- `fit_ram_twolevel`: v0.2 within/between moment decomposition
- `fit_ram_twolevel_ml`: exact complete-data Gaussian clustered likelihood
- `fit_ram_twolevel_fiml`: exact observed-subvector likelihood for arbitrary IEEE-NaN patterns, with an optional unrestricted within/between H1 likelihood

v0.6 also adds `fit_random_coefficient_ml`, `random_coefficient_loglik`, and `random_effects_eb`. This is a general Gaussian random-coefficient likelihood with multivariate outcomes, fixed covariates, arbitrary random intercept/slope design columns, unrestricted random-effect covariance, unrestricted residual covariance, and empirical-Bayes cluster random effects.

The complete-data exact RAM likelihood uses the determinant identity

`|I_m x Sigma_W + J_m x Sigma_B| = |Sigma_W|^(m-1) |Sigma_W + m Sigma_B|`

for each cluster, so latent structure at both levels can be supplied through separate RAM models without constructing giant cluster covariance matrices.

## Robust and diagnostic inference

- casewise and cluster-robust sandwich SEs
- covariance-structure Satorra-Bentler scaling
- `covariance_scaled_tests`
- `scaled_tests_from_ugamma`
  - mean-scaled Satorra-Bentler statistic
  - mean/variance-adjusted Satterthwaite statistic and adjusted df
  - scaled-shifted statistic and shift parameter
- `yuan_bentler_from_traces` for Yuan-Bentler scaling and H0/H1 scaling-factor bookkeeping
- `hayakawa_trace_corrected` and `hayakawa_adjusted_tests`, implementing the current lavaan/Srivastava-Himeno-Yamada-Hayakawa unbiased second-trace estimator from casewise moment vectors
- `browne_residual_test` and `browne_residual_nt` for Browne residual-based ADF/normal-theory chi-square calculations
- nuisance-adjusted modification indices/EPCs/joint score tests
- linear Wald tests
- nonparametric bootstrap inference

### v0.7 mixed Muthen / conditional-x path

`muthen1984_mixed` accepts numeric and ordered endogenous variables together with optional exogenous covariates. Numeric variables use conditional linear regressions; ordered variables use conditional probit regressions. Pairwise stage-2 correlations are estimated as residual Pearson, conditional polyserial, or conditional polychoric correlations according to variable type. The returned covariance metric rescales numeric variables by their conditional residual standard deviations. When `compute_gamma=.true.`, the complete statistic vector (univariate parameters plus correlations) gets a deterministic delete-one stagewise Gamma and WLS inverse. This is a practical parity path for mixed/conditional-x data; the all-ordinal `muthen1984_ordinal` routine remains the direct analytic A11/A21/A22 implementation.

### v0.7 missing random coefficients and robust differences

`random_coefficient_loglik_missing`, `fit_random_coefficient_missing_ml`, and `random_effects_eb_missing` evaluate only observed response elements within each cluster while retaining the exact random-intercept/slope covariance induced by the cluster design. `satorra_bentler_difference_2001` and `satorra_bentler_difference_2010` expose lavaan-style scaled nested-model difference calculations directly.

## Fit diagnostics and prediction

`sem_fit_result` includes parameter estimates, SEs, vcov, implied moments, objective, log likelihood, chi-square/df where available, AIC/BIC, RMSEA, CFI, TLI, and SRMR.

Also included:

- covariance residuals
- standardized RAM coefficients/covariances
- regression factor scores
- RAM simulation
- `vec`, `vech`, inverse `vech`
- duplication/commutation matrices
- symmetric square roots
- orthogonal complements
- SPD/general solves

## Build

```text
fpm build
fpm test
```

Dependencies are vendored as FPM path dependencies. No R installation is required.

## Current parity boundary

v0.7 covers the main standalone numerical workflows in EFA, SAM, MIIV, multilevel Gaussian SEM, categorical WLS/PML, MML, multigroup fitting, constraints, bootstrap, and robust inference. It now also includes mixed continuous/ordinal conditional-x statistics, missing-response random-coefficient likelihoods, automatic/name-aware MIIV marker rewriting, SAM joint block covariance and second-order bias propagation, and explicit Satorra-Bentler 2001/2010 difference-test formulas. Remaining differences are narrower: the exact analytic mixed/conditional-x Muthen NACOV (v0.7 uses stagewise delete-one Gamma for this branch), complete lavaan parameter-table constraint rewriting around MIIV, fully structural adaptive high-dimensional MML, estimator-specific robust-difference variants beyond SB 2001/2010, and R syntax/S4 infrastructure.

See `PORTING_NOTES.md` for exact-versus-narrower implementation details.
