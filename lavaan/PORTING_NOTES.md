# Porting notes

## Scope

Upstream lavaan 0.7-2 contains more than 132,000 lines of R. The Fortran project targets numerical SEM algorithms and matrix representations, not a line-by-line recreation of R syntax/S4 infrastructure.

## EFA

`efa_principal_axis` iterates squared multiple-correlation communalities to a principal-axis solution. `efa_ml` uses the standard concentrated Gaussian factor-analysis likelihood: for a proposed uniqueness diagonal, the whitened covariance eigenstructure supplies the concentrated loading solution. The uniquenesses are optimized with the library's BFGS routine.

Rotation is delegated to the vendored `GPArotation-fortran` translation. This is a direct numerical dependency rather than an approximate local replacement and provides the major orthogonal and oblique criteria used by lavaan.

The current EFA interface accepts a covariance matrix and factor count. lavaan's EFA block syntax, block-specific constraints, automatic starts, and S4 summary tables are not reproduced.

## SAM

`sam_fit_cov`/`sam_fit_data` implement the central local two-stage idea: estimate the measurement parameters, insert them as fixed values into a structural RAM template, then fit only the stage-2 structural map. This is useful and source-transparent, but it does not yet implement lavaan's complete local/global SAM variance propagation, bias correction, or block-combination bookkeeping. Stage-2 Hessian SEs therefore condition on the stage-1 measurement estimates.

## MIIV / 2SLS

`miiv_2sls` is equation-level raw-data 2SLS. It returns homoskedastic 2SLS covariance, first-stage F diagnostics, and a Sargan statistic. `miiv_2sls_cov` performs the same projected-moment solve from covariance blocks.

`ram_miiv_candidates` screens candidate RAM nodes using the population model condition `Cov(x_z, epsilon_y)=0`, computed as `[(I-A)^-1 S]_{z,y}`, and requires nonzero model-implied covariance with at least one endogenous predictor. This captures the numerical core of model-implied instrument screening, but it does not rebuild lavaan's full partable-driven equation system or categorical/equality-constrained MIIV machinery.

## Categorical Gamma matrix

v0.2's `ordinal_wls_correlation_weights` remains available for a correlation-only statistic vector. v0.3 adds `categorical_wls_statistics`, whose statistic vector contains all marginal thresholds followed by off-diagonal polychoric correlations.

The asymptotic covariance is estimated by delete-one jackknife and scaled to the usual `N * Var(s)` Gamma convention. This gives a coherent threshold/correlation WLS/DWLS matrix but is not a transcription of lavaan's analytic Muthen-style categorical Gamma recipes. The jackknife is substantially more expensive for large data sets.

## PML

`fit_ram_pml_ordinal` retains the all-ordinal pairwise likelihood from v0.2. `fit_ram_pml_mixed` extends the composite likelihood to continuous-continuous Gaussian pairs and continuous-ordinal conditional-normal pairs, while keeping empirical ordinal thresholds fixed.

This covers the main pair likelihood types, but not every lavaan PML test correction, missing-data pattern, or covariate option.

## MML

`fit_mml_ordinal_factor` implements a normal-ogive ordinal factor marginal likelihood with tensor Gauss-Hermite quadrature. It supports one to three standardized latent factors and an explicit free loading mask. Marginal thresholds are empirical probit thresholds and residual standard deviations are chosen so each latent response has unit variance.

This is the key numerical integration path behind categorical factor MML, but it is deliberately narrower than lavaan's general MML representation: arbitrary structural latent regressions, mixed numeric/ordinal outcomes, dummy latent variables, conditional covariates, and adaptive/high-dimensional quadrature are not yet included.

## Two-level fitting

`fit_ram_twolevel` remains the v0.2 moment decomposition. `fit_ram_twolevel_ml` is a new exact complete-data Gaussian cluster likelihood. If cluster `g` has size `m`, within covariance `Sigma_W`, and between covariance `Sigma_B`, it evaluates the likelihood using

`|Sigma_cluster| = |Sigma_W|^(m-1) |Sigma_W + m Sigma_B|`

and the corresponding within-deviation plus cluster-mean quadratic form. This supports separate RAM latent structures at both levels without forming an `(m*p) x (m*p)` covariance matrix.

The exact routine currently assumes complete continuous data and a shared overall mean from the between-level RAM model. lavaan's multilevel missing-data likelihood, saturated H1 construction, random slopes, and cluster-specific covariates remain outside this release.

## Robust test variants

`robust_ml_inference` still provides sandwich covariance and the covariance-structure mean-scaled Satorra-Bentler statistic. v0.3 adds `covariance_scaled_tests`, which forms `U Gamma` and reports both the mean-scaled statistic and the Satterthwaite mean/variance-adjusted statistic:

- `c_SB = tr(U Gamma) / df`
- `c_MV = tr((U Gamma)^2) / tr(U Gamma)`
- `df_MV = tr(U Gamma)^2 / tr((U Gamma)^2)`

This closes the principal mean/variance adjustment but does not enumerate every lavaan estimator-specific scaled/shifted/Yuan-Bentler/FML correction.

## Multigroup and constraints

`ram_group_spec%link` maps group-local free parameters to global parameter IDs. Reusing IDs imposes equality constraints. `fit_ram_multigroup_cov` and `fit_ram_multigroup_fiml` support covariance/raw-data estimation, unequal group sizes, and missing-data patterns.

`fit_ram_cov_constrained` uses repeated quadratic-penalty BFGS solves with a tangent-space covariance projection for equality constraints.

## Numerical dependencies

- `numDeriv-fortran`: Hessians and numerical derivative support.
- `pbivnorm-fortran`: bivariate normal rectangle probabilities.
- `GPArotation-fortran`: gradient-projection factor rotation.

All are vendored FPM path dependencies.

## Licensing

Upstream lavaan declares `GPL (>= 2)`. The vendored dependencies are compatible with GPL-2.0-or-later redistribution, so this combined package remains GPL-2.0-or-later.

## v0.4 refinement notes

### SAM propagated covariance

`sam_propagate_uncertainty` retains the stage-2 conditional Hessian covariance and additionally computes a local sensitivity matrix `J = d(theta_stage2)/d(theta_stage1)` by refitting stage 2 at centered perturbations of every measurement parameter. The propagated covariance is

`V_stage2,prop = V_stage2,conditional + J V_stage1 J'`.

This closes the main local two-stage uncertainty omission from v0.3. It is not lavaan's full global/block SAM bias-correction bookkeeping, and in orthogonal parameterizations `J` may legitimately be zero.

### MIIV equation discovery

`ram_miiv_equations` scans nonzero RAM regressions, forms one equation per dependent node, screens observed candidate instruments by the model-implied disturbance-exogeneity condition, excludes descendants of the dependent variable, requires relevance for an endogenous predictor, and reports equation-level identification. This is substantially closer to automatic MIIV construction than the v0.3 single-equation screener, but it deliberately uses RAM node indices instead of reconstructing lavaan's full named partable/marker-variable rewriting.

### Categorical Gamma

`categorical_wls_statistics_analytic` removes the O(N) leave-one-out refitting cost. Threshold equations use the exact probit influence equation `I(Y<=c)-Phi(tau)`, while each polychoric component uses its pair-likelihood score. The sensitivity matrix of the stacked estimating equations is differentiated numerically, and the asymptotic Gamma matrix is the sandwich `A^{-1} B A^{-T}`. Thus threshold/correlation cross-uncertainty is propagated jointly.

This is an influence-function implementation of the same asymptotic object, not a literal transcription of every closed-form Muthen NACOV expression in lavaan. The original jackknife path remains available as an independent check.

### Missing-data two-level likelihood and H1

`fit_ram_twolevel_fiml` forms each cluster covariance as `I_m kron Sigma_W + J_m kron Sigma_B`, selects only observed scalar entries, and evaluates their exact Gaussian subvector likelihood. This handles arbitrary IEEE-NaN patterns within and across variables/clusters.

When requested, an unrestricted H1 model is fitted with separate Cholesky-parameterized within and between covariance matrices plus an unrestricted mean vector. The result reports H1 log likelihood, likelihood-ratio chi-square, and the H1-H0 parameter-count df. The older determinant-identity complete-data routine remains faster when there is no missingness.

### Generalized MML

`fit_mml_mixed_factor` extends the ordinal MML path to mixed continuous/ordinal data, missing responses, arbitrary fixed latent means/covariances, and one to three latent dimensions. Continuous intercepts and residual scales are estimated jointly with the requested loading pattern. Ordinal thresholds remain empirical marginal probit thresholds.

Latent structural means/covariances can therefore be represented in the integration distribution, but they are currently supplied rather than jointly parameterized as a full RAM structural model. Adaptive/high-dimensional quadrature remains outside this release.

### Additional robust tests

`scaled_tests_from_ugamma` now also reports the scaled-shifted T3 form using `a=sqrt(df/tr((UGamma)^2))`, scaling `1/a`, and shift `df-a*tr(UGamma)`. `yuan_bentler_from_traces` records the Yuan-Bentler scaling plus separate H1/H0 trace scaling factors when supplied. Hayakawa casewise corrected traces and every lavaan estimator-specific wrapper remain future refinements.


## v0.5 refinement notes

### Hayakawa corrected robust traces

`hayakawa_trace_corrected` is a direct numerical translation of the current lavaan implementation of the Srivastava (2005) / Himeno-Yamada (2014) unbiased estimator used by Hayakawa (2018). For centered casewise saturated moment vectors `Zc`, it forms only the Gram matrix `M = Zc U Zc'`, so no square root of `U` is needed. It reports both `tr(U Gamma)` with the `N-1` covariance divisor and the corrected estimator of `tr((U Gamma)^2)`. `hayakawa_adjusted_tests` applies lavaan's corrected mean/variance-adjusted and corrected scaled-shifted formulas.

### MIIV covariance methods

The new `lavaan_miiv_variance` module ports the covariance-statistic linearization formulas used upstream for ULS, GLS, 2RLS and RLS MIIV estimation. In particular, the GLS/2RLS/RLS Jacobians include the derivative of the covariance-dependent weight matrix, not just the fixed-weight projection. The v0.5 regression test compares these Jacobians with centered finite differences of the corresponding refitted estimators on a nonsaturated covariance model. This is still a matrix-native API; lavaan's partable/name/marker rewriting remains outside the Fortran layer.

### Yuan-Chan SAM global correction

`sam_yuan_chan_test` exposes the matrix-level core of lavaan's current global SAM correction for a single assembled statistic block. It residualizes the saturated-statistic covariance as `D = I - Delta_gamma P`, `Gamma_tilde = D Gamma D'`, then applies the structural-only Satorra-Bentler projection using `Delta_theta` and `W`. Multigroup callers can assemble the corresponding block-diagonal matrices before calling this routine. Full lavaan block bookkeeping and global SAM bias corrections are still not reproduced.

### High-dimensional MML by deterministic QMC

`mml_mixed_loglik_qmc` and `fit_mml_mixed_factor_qmc` use deterministic Halton points transformed to standard normals and then through the supplied latent Cholesky factor. This retains the same mixed continuous/ordinal conditional likelihood as `mml_mixed_loglik` but avoids tensor-product Gauss-Hermite growth. The QMC path supports up to eight latent dimensions. It is a practical deterministic high-dimensional integration alternative, not adaptive Gauss-Hermite quadrature and not a transcription of every lavaan MML optimizer/integration strategy.


## v0.6 refinement notes

### Stagewise Muthen categorical NACOV

`muthen1984_ordinal` is a direct numerical implementation of the all-ordinal, no-covariate branch of lavaan's current Muthen (1984) three-stage machinery. Stage 1 uses the exact ordinal marginal log-likelihood threshold scores. Stage 2 uses analytic polychoric rectangle scores: Plackett's derivative for rho and analytic boundary derivatives for thresholds. Stage 3 constructs the same block-triangular bread layout `B = [[A11,0],[A21,A22]]`, the stacked case-score meat, and `B^{-1} INNER B^{-T}`. The returned `Gamma` is on lavaan's `N * acov(s)` scale.

This is materially closer to upstream than the v0.4 generic influence-equation routine and avoids numerical differentiation of the full stacked sensitivity matrix. The current direct path is all-ordinal and has no sampling weights, exogenous regressors, empty-category pseudo-thresholds, clustering, or mixed continuous/polyserial branch; those remain with the more general existing paths or as future extensions.

### Marker-variable MIIV rewriting

`ram_miiv_marker_equations` accepts explicit latent-to-marker node mappings. A latent structural equation is rewritten on the observed marker scale, including loading-ratio coefficient transformations. Candidate observed instruments are then screened against the covariance of the full composite proxy disturbance, not merely the latent disturbance, and descendants of the structural outcome are excluded. This closes the main numerical marker-variable transformation gap while deliberately avoiding lavaan's R parameter-table/name manipulation layer.

### Adaptive MML

`mml_mixed_loglik_adaptive` locates each case's latent posterior mode and numerical Hessian, uses its inverse as a local proposal covariance, and evaluates a posterior-centered Gauss-Hermite importance quadrature. It supports mixed continuous/ordinal data, missing values, arbitrary fixed latent means/covariances, and one to three latent dimensions. `fit_mml_mixed_factor_adaptive` exposes the same loading/intercept/residual-scale optimizer as the nonadaptive mixed MML path. The deterministic QMC path remains preferable for larger latent dimensions.

### Continuous SAM Gamma

`sam_continuous_gamma` provides the continuous ADF moment covariance for means and vech covariances. `sam_browne_unbiased_gamma` applies Browne's finite-sample unbiased covariance-moment correction, including the third-order mean/covariance rescaling used by lavaan. This supplies the finite-sample Gamma object needed by local/two-step SAM corrections without changing the existing stage-1/stage-2 fitting API.

### Random coefficients and slopes

`random_coefficient_loglik` evaluates the exact Gaussian cluster likelihood for multivariate outcomes with fixed design `X`, random design `Z`, unrestricted random-effect covariance `G`, and unrestricted residual covariance `R`. `fit_random_coefficient_ml` estimates fixed coefficients plus Cholesky-parameterized `G` and `R`. `random_effects_eb` returns posterior empirical-Bayes random effects and posterior covariance by cluster. This captures the central numerical random-intercept/random-slope likelihood and prediction machinery, but it is a matrix-native interface rather than lavaan's `rv()`/parameter-table mapping and currently assumes complete outcomes.

### Browne residual tests

`browne_residual_test` evaluates `N r' Gamma^{-1} r` for a supplied asymptotic covariance, and `browne_residual_nt` constructs the Gaussian normal-theory vech covariance from the model covariance before evaluating the residual statistic.


## v0.7 refinements

### Mixed continuous/ordinal Muthen path

The upstream mixed/conditional-x branch uses the full Muthen (1984) three-stage analytic derivative machinery. v0.7 now matches the stage-1 and stage-2 model definitions (linear/probit univariate fits and numeric-numeric/polyserial/polychoric latent correlations), including the covariance-metric rescaling. For the mixed branch Gamma, v0.7 deliberately uses deterministic delete-one stagewise refitting instead of reproducing every closed-form A21/H term. The all-ordinal path in `lavaan_muthen1984` remains analytic.

### Random coefficients with missing responses

The v0.7 missing-data random-coefficient likelihood constructs each cluster covariance only for observed response elements. Random-slope covariance contributions are unchanged; residual covariance contributes only when two observed elements come from the same row. This is the exact observed-subvector Gaussian likelihood under the fitted random-coefficient model.

### MIIV names and markers

Fortran does not reproduce lavaan's R parameter-table object. v0.7 adds automatic marker selection, a compact parameter-table marker descriptor, and named equation/instrument output while retaining the RAM-based disturbance screening introduced earlier. Equality/nonlinear parameter-table rewriting remains a higher-level frontend concern.

### Robust difference tests

The SB-2001 scaling factor is `(df0*c0 - df1*c1)/(df0-df1)`. The SB-2010 routine accepts the cross-model M10 scaling factor explicitly, matching lavaan's `lav_test_diff_sb2010` numerical formula without attempting to construct the R model object internally.
