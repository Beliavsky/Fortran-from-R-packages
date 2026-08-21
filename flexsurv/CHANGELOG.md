# Changelog

## 0.3.0

Numerical-parity release for flexsurv 2.3.2.

- Added predictable time-dependent covariate updates to clock-reset multi-state
  simulation. Parametric transition design columns can advance with absolute time;
  spline gamma coefficients can receive explicit predictable time slopes.
- Added a `splines2ns`-compatible natural cubic B-spline basis and derivative,
  including linear boundary extrapolation, and integrated it into spline fitting
  and gamma-interaction models.
- Added source-style initialization for the built-in parametric families. The
  Weibull, Weibull-PH, exponential, log-normal and log-logistic paths can use the
  supplied survival AFT (`survreg`) translation; the other families use the
  corresponding flexsurv moment formulas.
- Added selectable BFGS and Nelder-Mead optimization and a BFGS-to-simplex fallback.
- Reworked Louis observed information so multinomial membership scores/Hessians
  are analytic while component survival derivatives remain local numerical
  derivatives. Retained the v0.2 all-numerical routine as a regression reference.
- Added `flexsurv_shared_multistate` for one shared regression coefficient vector,
  cross-transition cumulative-hazard covariance and joint parameter bootstrap.
- Added five focused parity test programs. The strict suite now has 23 tests.
- Compiles the supplied `survival_aft` modules used for source-style starts.

## 0.2.0

Parity expansion for flexsurv 2.3.2.

- Added full mixture regression with multinomial-logit membership covariates,
  known/unknown/partially known event types, arbitrary component regressions,
  direct and EM fitting, and direct/Louis covariance estimates.
- Added `fmixmsm` pathway/final probabilities, means and simulated quantiles,
  including full-parameter bootstrap confidence intervals.
- Added rate-table-aware `standsurv` all-cause/attributable estimands,
  transformations, contrasts, delta-method and normal-bootstrap uncertainty.
- Added covariate effects on individual Royston-Parmar gamma coefficients.
- Added semi-Markov clock-reset simulation, `pmatrix.simfs`, LOS calculations,
  and parameter-bootstrap uncertainty.
- Added `pfinal_fmsm` and `simfinal_fmsm` final-state summaries and intervals.
- Added nonparametric Aalen-Johansen comparison helpers.
- Extended `flexsurvrtrunc` with arbitrary fixed distribution parameters and
  bootstrap-ready covariance; fixed coordinates remain exact in bootstrap draws.
- Integrated the translated relsurv rate-table engine into the compiled source.
- Expanded the strict regression suite from 10 to 18 test programs.

## 0.1.0

Initial modern Fortran/FPM computational translation of flexsurv 2.3.2.
