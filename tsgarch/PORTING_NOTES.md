# Porting notes

## Design

The R package uses S3 method dispatch, TMB automatic differentiation, nloptr,
R date-series classes, and parallel helpers. The Fortran port separates those
layers into explicit numerical types and procedures. All state needed to repeat
a calculation is stored in a `garch_spec`, `garch_parameters`, or result object.

## Model recursions

The translated filter implements the eight upstream volatility-model families:
GARCH, GJR-GARCH, APARCH, EGARCH, family GARCH, component GARCH, IGARCH, and
EWMA. Conditional means may be zero or constant. Variance regressors enter the
variance intercept, and variance targeting replaces the freely estimated
intercept with the sample-targeted effective intercept.

The implementation supports arbitrary ARCH/GARCH orders where meaningful.
Component-GARCH state vectors expose permanent and transitory components in
filter and simulation results.

## Parameter constraints

Natural parameters are packed for optimization and unpacked into typed model
parameters. Bounds depend on the model and innovation distribution.

- Nonnegative ARCH/GARCH coefficients are bounded explicitly.
- IGARCH uses a reduced parameterization and derives the final coefficient so
  that the integrated equality holds exactly.
- EWMA uses one decay coefficient and derives the shock coefficient.
- Stationarity, pairwise family-GARCH restrictions, and component-GARCH
  restrictions receive large objective penalties when violated.
- Distribution shape, skewness, and GH parameters use the domains from the
  vendored `tsdistributions` port.

## Optimization and derivatives

Upstream estimation combines TMB automatic differentiation and nloptr. This port
uses a bounded Nelder-Mead optimizer from the GPL-2-compatible
`tsdistributions` translation. Numerical central differences provide Hessians
and observation scores.

This choice has several consequences:

- Estimates should be numerically close when the same local optimum is found,
  but optimizer paths and termination messages are not expected to match R.
- Numerical standard errors will not be bit-for-bit equal to TMB automatic-
  differentiation results.
- The fit can be more sensitive to starting values for difficult high-order,
  heavy-tailed, or near-integrated specifications.
- `fit_options` exposes iteration, tolerance, simplex-scale, and inference
  controls instead of nloptr's algorithm-specific option list.

The user-supplied LGPL-3.0-or-later `nloptr` translation is not compiled or
linked because `tsgarch` is GPL-2.0-only. It is retained unchanged only as a
provenance attachment.

## Forecasts and simulation

`forecast_garch` uses conditional Monte Carlo simulation for every model. This
provides one common path for nonlinear recursions and non-Gaussian innovations.
Upstream code can use analytic or model-specific shortcuts in some cases, so
small Monte Carlo differences are expected. Increase `paths` when tail
quantiles or long horizons require greater precision.

The port currently provides parametric simulation. R-side residual/bootstrap
prediction modes and parallel execution are not implemented.

Random-number streams are deterministic for a given Fortran seed, but differ
from R's generator and therefore are not bit-for-bit reproducible across
languages.

## Rolling backtests

`backtest_var` supports:

- an expanding estimation sample when `window` is absent;
- a rolling fixed-size sample when `window` is present;
- user-selected refit intervals;
- one-step VaR forecasts;
- Kupiec unconditional coverage;
- Christoffersen independence;
- Christoffersen conditional coverage.

Dates, time zones, `xts`, and `zoo` metadata remain the caller's responsibility.

## Distribution dependency

The complete numerical `tsdistributions` Fortran dependency is vendored to keep
the package self-contained. Its GH-family routines in turn contain GPL-2.0-or-
later `ghyp` translations. File-level SPDX identifiers preserve those terms.

## Omitted R infrastructure

The following items do not represent independent numerical algorithms and were
not translated:

- S3 classes, print/summary methods, and formula parsing;
- `data.table`, `xts`, `zoo`, and `lubridate` adapters;
- plots and `flextable` formatting;
- `future` parallel scheduling and progress reporting;
- R list conversion helpers such as `to_multi_estimate`;
- compiled TMB/Rcpp registration wrappers.
