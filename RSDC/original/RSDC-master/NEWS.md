# Changes in Version 1.7-0 (DA,BS,RN)
- Bug fix: `rsdc_viterbi()` no longer errors on a single-observation input
  (the forward recursion indexed out of bounds).
- `rsdc_minvar()`/`rsdc_maxdiv()` validate their inputs before reconstructing
  the covariance: volatilities must be finite and positive and correlations
  finite with modulus below one, so a malformed column now raises an
  informative error naming the offending row instead of silently producing a
  wrong-signed covariance matrix.
- `rsdc_maxdiv()` checks that `Rsolnp::solnp()` actually converged (and that
  its weights are finite and sum to one) rather than accepting whatever it
  returns. Both optimizers now count the periods in which they fell back to
  equal weights, report the count as `$n_fallback`, and warn once if it is
  positive; previously the fallback was silent.
- `rsdc_forecast_ahead()` validates `horizon`, and `rsdc_forecast()`
  validates `control$threshold`, instead of failing later with a cryptic
  message.
- Bug fix (standard errors): the numerical Hessian is now rejected when the
  optimum sits within the finite-difference step of the feasible boundary
  (a correlation within ~2e-3 of +/-1, or a nearly singular correlation
  matrix). The objectives return a large finite penalty there, so the
  differenced Hessian was huge but finite and its inverse produced absurd,
  silently reported standard errors (e.g. 2e-08). `vcov()`/`summary()`/
  `confint()` now return `NA`/`NULL` with an explanatory warning instead.
- `vcov(type = "opg")` no longer requires the Hessian-based covariance: OPG
  needs only the per-observation scores, which are now stored whenever
  `compute_se = TRUE`. Scores are still withheld (and all covariance
  estimators refused) when a `"noX"` transition row had to be projected onto
  the simplex, since that point is not a stationary optimum.
- `Depends` is now `R (>= 4.0)`: the package compiles C++11 (from R 4.0 the
  default C++ standard) and declaring `CXX_STD` explicitly would raise an
  R CMD check NOTE.
- `rsdc_estimate()` now rejects fewer than two observations instead of
  returning a degenerate fit.
- Windows build fix: `src/Makevars`/`src/Makevars.win` now link BLAS/LAPACK
  explicitly (`$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`), as required by
  RcppArmadillo. Without it, `R CMD INSTALL` on Windows failed at the link
  step with undefined references (`dgemv_`, `dgemm_`, `dsyev_`, ...); macOS
  linkers tolerate the undefined symbols, which is why the omission was
  invisible on the development machine.
- The package article (methodology, a strictly out-of-sample five-industry
  study, and a Monte Carlo validation) is developed separately and is
  reproducible from data shipped with the package; it is no longer shipped as a
  package vignette. See `citation("RSDC")` for the current reference. The
  package ships the lightweight "Getting started" vignette (`rsdc-intro.Rmd`).
- Reparameterized global search: `rsdc_estimate()`'s DEoptim stage now
  explores canonical partial correlations (Joe, 2006, J. Multivariate
  Analysis 97(10), 2177-2189) instead of raw pairwise correlations. Every
  point of the search box maps to a valid positive-definite correlation
  matrix, so the out-of-the-box global search is feasible at any K (the
  natural box has essentially no feasible point from K ~ 5 on). For
  `method = "noX"` with N >= 3 the transition head is searched as bounded
  per-row softmax logits (roughly the historical [0.01, 0.99] box on the log-odds
  scale), which removes the row-simplex penalty wall while keeping every
  regime reachable; the N = 2 head and all natural-parameter results are
  unchanged. Within each search, the top three mutually distant members of
  the final DE population are refined by L-BFGS-B and the best refined
  point is kept (the best DE point is not always in the best basin).
- Replication diagnostic: with `control$n_starts >= 2`, a warning is issued
  when the best log-likelihood is reached by only one of the independent
  searches (tolerance 2.5 log-points), signaling a possibly non-global
  optimum. `control = list(n_starts = 4, cores = 4)` is the recommended
  protocol for multimodal surfaces (e.g. `"noX"` with N >= 3); the default
  remains a single search.
- `rsdc_starts()` is unchanged and remains a fast deterministic complement:
  when both routes agree, the optimum is cross-certified.
- Behavior change: `cores > 1` on a plain single fit no longer enables the
  internal `parallelType` of DEoptim; parallelism lives across the
  independent searches of the multi-start forms.
- New dataset `ff5ind`: daily value-weighted returns of five Fama-French
  industry portfolios (Manuf, Enrgy, HiTec, Hlth, Utils; Fama & French,
  1997, J. Financial Economics 43(2), 153-193; data (c) Fama and French,
  Kenneth R. French Data Library) with the MCCC Aggregate index (2025
  Sentometrics release, monthly, forward-filled across trading days) and
  the CBOE VIX (via FRED, VIXCLS), aligned on 5155 U.S. trading days,
  2005-01-03 to 2025-06-30. Reproducible build in `data-raw/ff5ind.R`.

# Changes in Version 1.6-0 (DA,BS,RN)
- New `rsdc_starts()`: data-driven warm starts for high-dimensional problems.
  In moderate dimensions the DEoptim global search works well, but for K >= 5
  a random vector of pairwise correlations is almost never positive definite,
  so its initial population contains no feasible point. `rsdc_starts()` builds
  feasible starting vectors from the empirical correlations of low/…/high
  average-correlation sub-periods (positive definite by construction), with
  `n_starts` shrinkage variants as protection against local optima.
- `rsdc_estimate()` accepts an `rsdc_starts` object via `control$start`: each
  start is refined by the local optimizer (the global search is skipped), the
  highest-likelihood fit is kept, and the returned object carries
  `start_logliks` (spread across starts) and `start_pars` (per-start optima).
  `control$cores` parallelizes both multi-start forms (`n_starts` seeds or
  `rsdc_starts` warm starts); for a plain single fit it is forwarded to the
  backends, where `cores > 1` enables parallel DEoptim (`parallelType = 1`),
  as before.
- New `cores` argument for `rsdc_bootstrap()` and `rsdc_corr_bands()`:
  task-level parallelism via the base `parallel` package (forked workers on
  Unix-alikes, a PSOCK cluster on Windows). All random draws are generated up
  front from one uninterrupted RNG stream and each parallelized task is
  deterministic given its inputs, so results are identical for any number of
  cores; only the wall time changes.
- Bug fix: the estimation backends no longer call `set.seed()` when a warm
  start skips the global search. Previously every warm-started refit reset the
  global RNG stream, so `rsdc_bootstrap()` — which interleaves simulation and
  warm-started re-estimation — produced nearly identical replicates from the
  second one on, drastically understating bootstrap standard errors and
  intervals (also via `vcov(type = "bootstrap")` and
  `confint(type = "bootstrap")`).

# Changes in Version 1.5-0 (DA,BS,RN)
- Bug fix (audit): `rsdc_forecast()` now returns a scalar `NA` BIC (instead of a
  silent zero-length `numeric(0)`) when `final_params$log_likelihood` is missing.
- Bug fix (audit): a warm start (`control$start`) lying outside the default
  L-BFGS-B search box — as produced by the identifiability relabeling for
  `N >= 3` (re-referenced softmax betas, former row-complement transition
  entries) — is no longer silently projected onto the box; the box is widened
  element-wise to keep the supplied start feasible (affects multi-start refits
  and `rsdc_bootstrap()` warm starts).
- CRAN hygiene: stray local session files under `inst/simulation/results/` are
  now build-ignored, and the declared `Rdpack` import is used
  (`importFrom(Rdpack, reprompt)`), clearing both `R CMD check` NOTEs.
- `rsdc_minvar()`/`rsdc_maxdiv()`: when `value_cols` is a character vector and
  `y` carries column names, misaligned columns are now detected — reordered to
  `value_cols` (with a warning) when the names match as a set, otherwise a
  warning states that positional order is assumed. Previously the names were
  silently overwritten.
- `rsdc_forecast_ahead()` now validates that `X_future` has the same number of
  columns as the estimation covariate matrix `X`.
- `plot.rsdc_fit()` error message corrected: stored regime probabilities do not
  depend on `compute_se`.
- Documentation: `rsdc_corr_bands()` documents the truncation of invalid
  Gaussian draws near parameter bounds.
- Removed dead `requireNamespace("mvtnorm")` guards (`mvtnorm` is a hard
  dependency in `Imports`).
- Multi-start estimation: `control$n_starts` repeats the global+local search from
  several seeds and keeps the highest-likelihood fit, storing `start_logliks`.
  Since `DEoptim` is a stochastic global search, this is primarily a stability
  diagnostic across seeds (is the optimum reproducible?), not a replacement for it.
- `rsdc_forecast_ahead()`: genuine multi-step-ahead forecasts of the regime
  distribution and implied correlations, propagating the terminal filtered state
  through the Markov chain (`X_future` for the time-varying case).
- `rsdc_corr_bands()`: pointwise uncertainty bands for the predicted correlation
  path, by drawing parameters from the asymptotic sampling distribution and
  re-running the filter.
- `rsdc_viterbi()`: most likely (MAP) regime path via the Viterbi algorithm.
- broom tidiers `tidy()`, `glance()`, `augment()` and a `ggplot2::autoplot()`
  method for `rsdc_fit` (ggplot2 is a soft dependency in `Suggests`).
- Documentation: the "Getting started" vignette gains a section demonstrating the
  new features above (multi-start, bands, Viterbi, multi-step forecasts, tidy
  output), and the broom/autoplot methods gain runnable examples.

# Changes in Version 1.4-0 (DA,BS,RN)
- Performance: the Hamilton-filter log-likelihood used during estimation is now
  evaluated in C++ (Rcpp/RcppArmadillo), matching the R reference to ~1e-8 and
  cutting estimation time by roughly an order of magnitude. The pure-R
  `rsdc_hamilton()` is retained (and used as the equivalence reference).
- Optimizer control: `control` now also accepts `cores` (parallel `DEoptim` via
  `parallelType = 1`) and `start` (a warm-start parameter vector that skips the
  global search and goes straight to local refinement).
- Parametric bootstrap: `rsdc_bootstrap()` simulates from the fitted model,
  re-estimates each replicate warm-started at the MLE, and returns bootstrap
  standard errors / percentile intervals that respect the parameter bounds.
- Robust standard errors: `vcov()`/`confint()` gain a `type` argument —
  `"hessian"` (default), `"opg"` (outer product of gradients), or `"sandwich"`
  (QML). Per-observation scores are stored on the fitted object.
- `summary()` now reports delta-method standard errors on the natural scale of
  the transition probabilities, plus regime diagnostics (stay probability,
  expected duration, ergodic distribution).
- Added a `plot()` method for fitted models (smoothed/filtered regime
  probabilities).
- Fitted objects now store the estimation residuals/covariates and the
  filtered/smoothed regime probabilities.
- Input-validation hardening (audit): `rsdc_estimate()`/`rsdc_hamilton()` now
  reject single-series input (`K < 2`) and a covariate matrix `X` whose row
  count does not match the data; `rsdc_hamilton()` validates that a user-supplied
  fixed transition matrix `P` is row-stochastic (finite, non-negative, rows
  summing to 1).
- `confint()` now uses the same safe square-root as `summary()`, so non-positive
  variance estimates map to `NA` instead of producing inconsistent `NaN`
  intervals.
- For the fixed-transition `noX` model with `N >= 3`, a transition row that the
  optimizer leaves infeasible (free probabilities summing above 1) is still
  projected onto the simplex, but the fit is now flagged as non-converged and its
  standard errors are suppressed (the projected point is not a valid basis for
  inference).

# Changes in Version 1.3-0 (DA,BS,RN)
- `rsdc_estimate()` now returns an object of class `"rsdc_fit"` with standard S3
  methods: `print`, `summary`, `coef`, `logLik`, `nobs`, `vcov`, `confint`,
  `predict`, and `simulate`. `AIC()`/`BIC()` work out of the box via `logLik`.
- Standard errors: the estimator now returns the observed-information
  variance-covariance (`vcov`) and standard errors (`summary`/`confint`), computed
  from the numerical Hessian of the negative log-likelihood at the MLE.
- Arbitrary number of regimes: `N >= 4` is now supported for `"noX"` and `"tvtp"`
  (previously capped at `N = 3`).
- `control` now forwards optimizer settings (`itermax`, `NP`, `parallelType`,
  `steptol`, `maxit`, and `compute_se`) to `DEoptim`/`optim`, enabling faster runs.
- Robustness: warns when a transition probability is pinned to its bound and when
  a `"tvtp"` covariate matrix `X` has no intercept/constant column.
- Added a "Getting started" vignette.
- Moved the 596 KB source workbook out of `inst/extdata` into build-ignored
  `data-raw/` to shrink the package tarball.

# Changes in Version 1.2-0 (BS,DA,RN)
- Added `mccc` dataset: daily Media Climate Change Concerns (Aggregate) index,
  forward-filled from monthly and aligned row-for-row to `greenbrown`. Supplies
  the exogenous TVTP covariate used in the vignette, so the empirical workflow
  is fully reproducible from packaged data alone (no external Excel files).
- Added lightweight raw source `inst/extdata/mccc-monthly.csv` and builder
  `data-raw/mccc.R`.
- Fixed `greenbrown` documentation: date range is 2014-01-02 to 2022-12-30.

# Changes in Version 1.1-2 (BS,DA)
- DESCRIPTION improved following the CRAN guidelines
- Seed is now a control parameter
- Removed examples for unexported functions
- Fixed do-not-run examples
- Added green-brown-ptf data in extdata
- Data added under the right .rda format
- Data Documentation added

# Changes in Version 1.1-1 (DA)
- URLs fixed

# Changes in Version 1.1-0 (DA)
- First release public version
