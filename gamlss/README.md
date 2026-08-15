# gamlss-fortran v0.9.0

Modern Fortran 2018/FPM translation of the computational core of the R package
`gamlss` 5.5-0 (Rigby, Stasinopoulos and contributors).

This is a matrix-first numerical library. R formula/model-frame/S3/S4 machinery
is replaced by explicit design matrices, offsets, weights and penalty matrices.
Plotting code is intentionally omitted.

## Core fitting engine

- Rigby-Stasinopoulos (RS) parameter-wise Fisher/IRLS iteration.
- Cole-Green (CG) joint iteration with cross-parameter observed-information
  blocks.
- Mixed RS -> CG fitting.
- Separate design matrices for `mu`, `sigma`, `nu`, and `tau`.
- Observation weights and parameter offsets.
- Fixed parameters, parameter-specific step lengths and automatic step halving.
- Quadratic penalties and parameter-specific smoothing parameters.
- Optional ML-style in-iteration estimation of a scalar penalty multiplier for
  each distribution parameter.
- Effective degrees of freedom, penalized deviance, AIC/GAIC and SBC/BIC.
- Prediction of all fitted distribution parameters on new design matrices.

## New in v0.9.0

### Upstream `getMarginal()` random-effect prediction

`marginal_predict_eta` and `get_marginal_random_intercept` translate the
computational behavior of upstream `getMarginal()` in `R/random.R`.  The four
methods are exposed as `MARGINAL_INTEGRATE`, `MARGINAL_QFUNCTION`,
`MARGINAL_RANDOM`, and `MARGINAL_NONE`.

- `MARGINAL_INTEGRATE` numerically integrates the parameter inverse link over a
  zero-mean Gaussian random effect with fitted standard deviation `sigma_b`.
- `MARGINAL_QFUNCTION` uses exactly the upstream deterministic grid
  `qnorm(0.001, 0.002, ..., 0.999) * sigma_b` and averages the inverse link over
  all 999 points.
- `MARGINAL_RANDOM` averages independent Gaussian random-effect draws; its
  default is the upstream value of 10,000 draws per prediction.
- `MARGINAL_NONE` removes the fitted random-term contribution and applies the
  inverse link without marginalization.

The low-level routine accepts the random-effect-free linear predictor directly,
which is the matrix-first analogue of upstream `etam <- lp - rt`.  The object
adapter accepts `random_intercept_result_t` plus the group vector and performs
that subtraction automatically.  The random-intercept fit now stores
`random_edf` and `sigma_b`; `sigma_b` follows the upstream `random()` estimate
`sqrt(sum(b**2) / edf)`.

The implementation uses the fixed family links in the translated
`gamlss.dist` backend.  Upstream R can construct families with alternative
runtime link choices; that dynamic family-object feature remains outside the
matrix-first Fortran API.

## New in v0.8.0

### Adaptive importance-sampling marginal random-effects likelihood

`fit_gamlss_joint_random_effects_ais` adds a moderate-dimensional alternative
to both the scalable Laplace/posterior-moment fit and the low-dimensional tensor
Gauss-Hermite fit.  It first obtains a GAMLSS/Laplace starting fit, then for
each group finds the conditional posterior mode of the full joint random-effect
vector with a damped Newton solve.  The inverse local posterior Hessian defines
an adaptive Gaussian proposal.  Deterministic Halton quasi-Monte-Carlo normal
points with antithetic pairing are then importance weighted against the exact
conditional GAMLSS likelihood times the Gaussian random-effect density.

The QMC cost is approximately linear in `qmc_points` rather than exponential in
latent dimension.  The implementation accepts combined active random-effect
dimensions up to 16; this is a practical numerical guard, not a tensor-grid
limit.  The returned result includes group marginal log likelihoods, posterior
means, full posterior covariance matrices, effective sample sizes, and the
minimum ESS as an integration-quality diagnostic.

By default AIS refines the marginal likelihood/posterior moments at the
Laplace starting estimates.  `refine_parameters=.true.` additionally optimizes
the fixed effects and full random-effect covariance under the deterministic AIS
marginal likelihood.  The inner posterior-mode solver is separate from the
outer BFGS optimizer, so this path does not require recursive optimizer calls.

`fit_gamlss_joint_random_effects_ghq` remains preferable for very small latent
blocks when high-order tensor quadrature is affordable; the v0.6 Laplace path
remains preferable for very large or latency-sensitive models.

## New in v0.7.0

### Discrete and mixed Gaussian-copula likelihoods

`fit_gamlss_gaussian_copula_mixed` extends the v0.6 copula fitter to
discrete families and mixed atomic/continuous families such as `BEINF`.  For
each correlation group it partitions coordinates into continuous observations
and atomic observations.  Continuous coordinates are conditioned on exactly;
atomic coordinates contribute latent-normal rectangle probabilities with
bounds `Phi^-1(F(y-))` and `Phi^-1(F(y))`.  Therefore the implementation
reduces to the ordinary Gaussian-copula density when all coordinates are
continuous, and to a multivariate-normal rectangle probability when all are
discrete.

`mvn_rectangle_probability` evaluates the rectangle probability with the Genz
sequential-conditioning transformation and deterministic Halton/antithetic
quasi-Monte Carlo points.  `mvn_conditional` and `mvn_logpdf` expose the
underlying Gaussian conditioning and density kernels.  The number of QMC points
is user-selectable so likelihood accuracy can be traded against fitting cost.

### Gauss-Hermite marginal random-effects likelihood

`fit_gamlss_joint_random_effects_ghq` directly integrates the joint Gaussian
random-effect vector instead of profiling random effects and applying a
Laplace/posterior-moment covariance update.  It supports a full cross-parameter
covariance and arbitrary q-column random-effect designs when the combined
latent dimension is at most four.  Five- and seven-point Gauss-Hermite rules
are available per latent dimension.  The result retains the marginal log
likelihood, full covariance, fixed-parameter covariance, group log likelihoods,
and quadrature posterior means of the random effects.

The v0.6 `fit_gamlss_joint_random_effects` remains the preferred scalable path
for larger latent blocks.

## New in v0.6.0

### Gaussian-copula joint likelihood for continuous margins

`fit_gamlss_gaussian_copula` adds a genuine joint-likelihood alternative to
`fit_gamlss_correlated_rs`.  For continuous GAMLSS margins it maximizes

```text
sum_i log f_i(y_i) + log c_R(F_1(y_1), ..., F_n(y_n))
```

using a Gaussian copula and the translated `nlme` correlation structures.
Marginal regression coefficients and non-fixed correlation parameters are
optimized jointly.  Grouped time/space correlation blocks are supported.
The result retains marginal, copula and joint log likelihoods, Gaussianized
scores, correlation parameters and the joint parameter covariance.  This API
intentionally rejects discrete families and the mixed-mass `BEINF` family,
because their exact copula likelihood requires multivariate rectangle
probabilities rather than a continuous copula density.

### Cross-parameter random-effect covariance

`fit_gamlss_joint_random_effects` fits grouped random effects with one full
covariance matrix over every active `(distribution parameter x random term)`
combination.  For example, random intercept/slope terms on both `mu` and
`sigma` produce a 4 by 4 covariance containing within-parameter and
cross-parameter terms such as `Cov(b_mu0,b_sigma0)` and
`Cov(b_mu1,b_sigma1)`.  The coefficient step maximizes the joint penalized
GAMLSS likelihood; covariance updates use posterior second moments from the
optimizer curvature and are geometrically damped.

The older `fit_gamlss_multi_random_effects` remains available when independent
covariance blocks by distribution parameter are desired.

## New in v0.5.0

### Correlated non-Gaussian RS working responses

`fit_gamlss_correlated_rs` embeds the supplied `nlme` residual correlation
matrices in the parameter-wise RS/Fisher working-response updates.  The local
working covariance is

```text
D(1/sqrt(w_work)) V R V D(1/sqrt(w_work))
```

so IRLS curvature weights and an `nlme` correlation structure are both
retained.  This makes AR(1), CAR(1), ARMA, compound-symmetric, spatial and
unstructured correlations usable with translated non-Gaussian GAMLSS
families.  The shared correlation is estimated on the location working
response when it is not fixed.  A supplied variance-function specification is
evaluated as a base scale inside the working covariance; v0.5 does not jointly
re-estimate non-Gaussian variance-function parameters.

### Simultaneous random intercept/slope blocks

`fit_gamlss_multi_random_effects` allows arbitrary q-column grouped random
effect designs on several distribution parameters simultaneously.  Each
active parameter has its own q by q within-group covariance matrix, so a
correlated random intercept/slope pair can be fitted on both `mu` and `sigma`
in one model.  Covariance matrices are updated from conditional effects plus
posterior covariance blocks.  Cross-parameter random-effect covariance is not
imposed in this release.

### Validation and residual diagnostics

- `cross_validate_gamlss` performs explicit-fold out-of-sample validation and
  returns per-case/per-fold log likelihoods and fitted distribution parameters.
- `randomized_quantile_residuals_all` uses the generic family-support
  dispatcher and therefore covers the complete translated backend catalog,
  rather than only the early subset handled by the original v0.1 residual
  helper.

## New in v0.4.0

### Gaussian residual correlation and variance functions

`fit_gamlss_no_gls` provides an exact Gaussian/NO location-model adapter to
the supplied `nlme` GLS backend.  It supports the translated correlation
structures AR(1), continuous AR(1), ARMA, compound symmetry, exponential,
Gaussian, linear, ratio, spherical and unstructured correlation, together
with constant, fixed, group-specific, power, exponential, constant-plus-power
and constant-proportion variance functions.  The returned
`correlated_no_result_t` retains both the complete `gls_result` and a
GAMLSS-style result view.

### Random effects on multiple distribution parameters

`fit_gamlss_multi_random_intercept` can put independent grouped random
intercepts on any subset of `mu`, `sigma`, `nu`, and `tau` simultaneously.
Each active parameter receives its own quadratic penalty and independently
estimated RS variance ratio/penalty multiplier.  This complements the v0.3
`fit_gamlss_random_effects`, which supports correlated random slopes for one
target parameter at a time.

### Parameter-general stepwise GAIC

`stepwise_gaic_parameter` extends matrix-based forward/backward/bidirectional
selection to any available distribution parameter rather than only `mu`.
Candidate columns can therefore be selected for scale, shape, or kurtosis
equations while fixed design matrices are retained for the other parameters.

### Numerical diagnostics

- `worm_plot_diagnostics` computes equal-count local detrended-Q-Q cubic worm
  coefficients without plotting.
- `influence_from_hat` returns leverage-adjusted residuals and a Cook-style
  influence distance.
- `jarque_bera_statistic` provides a compact residual normality diagnostic.

### Executable-stack portability cleanup

Optimizer and integration callbacks in the GAMLSS censoring code and vendored
`gamlss.dist` numerical routines were refactored from nested procedures to
module procedures.  GNU Fortran no longer emits trampoline/executable-stack
linker warnings in the clean validation build.

## New in v0.3.0

### General grouped random effects

`fit_gamlss_random_effects` extends the v0.2 random-intercept wrapper to an
arbitrary grouped random-effect design such as random intercepts plus slopes.
The group-effect covariance may be diagonal or fully correlated.  The fitting
loop alternates between the ordinary GAMLSS RS/CG fit and an EM-like posterior
second-moment covariance update; for Gaussian `mu` models the supplied `nlme`
port can initialize the unstructured covariance matrix.

The returned `random_effects_result_t` contains the fixed/distributional GAMLSS
fit, group-by-random-term conditional effects, covariance and precision
matrices, compressed group levels, and outer covariance iterations.

### Delayed entry and survival adapters

`fit_gamlss_censored` now accepts an optional `entry` vector for left-truncated
(delayed-entry) likelihoods.  Exact, left-, right-, and interval-censored
contributions are conditioned on survival beyond the entry point for continuous
and discrete families.

New helpers are:

- `surv_interval2` for finite/infinite interval endpoints;
- `surv_counting_process` for `(start, stop, event)` data and delayed entry;
- `truncated_censored_case_loglik` for individual conditional likelihoods.

### Additive and multidimensional P-splines

- `build_additive_p_splines` combines multiple persistent P-spline terms into
  one identifiable design/penalty pair using sum-to-zero coefficient contrasts.
- `predict_additive_p_splines` reconstructs the identical basis at new points.
- `tensor_p_spline_2d` and `predict_tensor_p_spline_2d` provide row-wise tensor
  product bases with anisotropic Kronecker-sum penalties.

### Stepwise selection, bootstrap and intervals

- `stepwise_gaic_mu` supports forward, backward and bidirectional GAIC/BIC
  search from an explicit matrix scope.
- `bootstrap_gamlss_cases` performs nonparametric case-resampling refits for any
  translated family and retains parameter-coefficient samples and deviances.
- `bootstrap_percentile_ci` returns coefficient percentile intervals.
- `profile_likelihood_ci` converts a numerical coefficient profile into the
  usual one-degree-of-freedom likelihood-ratio confidence interval.

## New in v0.2.0

### Random effects

`fit_gamlss_random_intercept` integrates a grouped random intercept into any
selected GAMLSS parameter. The group coefficients are fitted as a quadratic
penalty block and the variance-ratio/penalty multiplier is updated inside the
RS iteration. For Gaussian location models the supplied `nlme` translation can
optionally provide the initial variance ratio.

The returned `random_intercept_result_t` contains the ordinary GAMLSS result,
compressed group levels, conditional group effects, fitted penalty multiplier,
and its working-scale variance ratio.

### Censored responses

`fit_gamlss_censored` implements exact, left-, right-, and interval-censored
likelihood contributions using the translated family PDF/CDF routines.
`family_cdf` covers the 62 generic families exposed by the vendored
`gamlss.dist-fortran v0.3.0` backend.

Censoring codes are:

```text
CENS_EXACT
CENS_LEFT
CENS_RIGHT
CENS_INTERVAL
```

`surv_right_censoring` converts the common `(time,event)` right-censored
representation used by `Surv(time,event)` into these arrays.

### Additional smoothers

- Persistent P-splines and natural splines from v0.1 remain available.
- `fractional_polynomial_basis`
- `select_fractional_polynomial`
- `predict_fractional_polynomial`
- 1D/2D local-polynomial `fit_loess` / `predict_loess`
- optional LOESS target-EDF span search
- `varying_coefficient_p_spline`
- `fit_monotone_p_spline` / `predict_monotone_p_spline`
- cyclic and ordinary finite-difference penalties
- ridge penalties

The fractional-polynomial search uses the upstream power grid
`(-2,-1,-0.5,0,0.5,1,2,3)` and repeated-power `x^p log(x)` rule.

### Penalized categorical fusion

`fit_pcat` translates the numerical core of upstream `pcat()`:

- all pairwise category differences;
- iterative adaptive Lp weights;
- fixed or ML-updated penalty multiplier;
- effective degrees of freedom;
- `pcat_fused_groups` for converting nearly equal fitted category effects into
  fused groups.

### Selection and profiling

- `forward_gaic_mu` performs forward GAIC/BIC column selection for the `mu`
  design while retaining the other GAMLSS parameter models.
- `profile_gamlss_coefficient` profiles any `mu`, `sigma`, `nu`, or `tau`
  regression coefficient over a supplied grid by fixing that coefficient via
  an offset and refitting the remaining model.

These are numerical matrix counterparts to the computational portion of the R
stepwise/profile workflows; formula/terms object generation remains outside the
Fortran layer.

## Existing v0.1 functionality

- P-spline design matrices with persistent knots and prediction.
- Natural-spline bases and cyclic penalties.
- `all_pair_difference_matrix`, translated from upstream `genD.c`.
- Penalized hat values / EDF.
- Randomized quantile residuals, residual moments and residual ACF.
- Deviance-increment comparison and likelihood-ratio statistic.
- Intercept-only distribution comparison (`fitDist`/`chooseDist` numerical role).
- LMS/centile fitting and prediction with BCCG, BCT, and BCPE.

## Distribution backend

The archive vendors `gamlss.dist-fortran v0.3.0`. It supplies 62 generic
fitted-family constants plus fixed-denominator DBI/ZIBB/ZABB APIs and their
distribution functions.

## Supplied dependencies

The supplied translations are retained as FPM path dependencies:

- `splines-fortran v0.1.0`
- `nlme` Fortran port 3.1.170
- `survival-fortran v0.1.0`

The GAMLSS code directly uses the spline backend and calls the supplied `nlme`
implementation both for optional Gaussian random-effect starts and for exact
Gaussian residual correlation/variance-function fits through `fit_gamlss_no_gls`.
The survival source is compiled as part of clean validation; GAMLSS censoring
itself uses a lightweight matrix representation rather than reproducing R's
`Surv` class.

## Validation

Run:

```text
./scripts/run_gfortran_tests.sh
```

The script compiles all declared dependency source trees and the GAMLSS library
with GNU Fortran 2018 and runtime bounds checking, then runs:

```text
test_core: PASS
test_lms_selection: PASS
test_v02: PASS
test_v03: PASS
test_v04: PASS
test_v05: PASS
test_v06: PASS
test_v07: PASS
test_v08: PASS
test_v09: PASS
```

The v0.2 suite tests censored-normal recovery, `nlme`-initialized grouped random
intercepts, fractional-polynomial power selection, LOESS reconstruction,
varying-coefficient and monotone P-splines, adaptive category fusion, forward
GAIC selection, and coefficient profile likelihood.

The v0.5 suite adds non-Gaussian correlated Gamma working-response fitting,
simultaneous `mu`/`sigma` random intercept/slope covariance recovery, K-fold
log-score validation, and full-family SHASH quantile residuals.

The v0.3 suite adds correlated random intercept/slope recovery, delayed-entry
likelihood recovery and survival adapters, additive/tensor P-splines,
backward/bidirectional selection, case bootstrap intervals, and profile-LR
confidence intervals.

The v0.4 suite adds AR(1) Gaussian residual-correlation recovery, simultaneous
`mu`/`sigma` random-intercept fitting, `sigma`-equation BIC selection, and
worm/influence/Jarque-Bera diagnostics.  Validation also verifies that the
previous GNU executable-stack warning has been removed.

FPM itself was not installed in the translation environment, but the manifest
uses the standard FPM path-dependency layout and the full declared source set is
compiled directly by the validation script.

## Scope still intentionally excluded

R formula/model-frame/terms/contrast parsing, S3/S4 object infrastructure,
printing, plotting, parallel stepwise wrappers, and graphical presentation
remain omitted.  Family links are fixed by the translated `gamlss.dist`
backend rather than selected dynamically through R family objects.  The
matrix GAIC selector does not recreate R formula scope/update objects.

Several later Fortran numerical extensions (copula likelihoods, GHQ and AIS
joint random-effect integration) go beyond the algorithms provided by upstream
`gamlss`; they are kept as explicitly named optional APIs.  v0.9's marginal
prediction routines, by contrast, are direct translations of upstream
`getMarginal()` behavior.

See `docs/TRANSLATION_NOTES.md` and `docs/API_MAP.md` for the exact mapping.
