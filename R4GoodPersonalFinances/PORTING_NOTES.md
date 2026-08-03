# Porting notes

## Numerical conventions

- Real calculations use `dp = kind(1.0d0)`.
- Gompertz truncation follows the upstream conditional-survival formula and
  returns zero at a degenerate maximum-age endpoint instead of NaN.
- Upper incomplete gamma values are evaluated by adaptive transformed
  quadrature, including negative shape parameters needed by the ruin formula.
- The gamma CDF in the ruin calculation uses standard series and continued
  fraction expansions.
- Portfolio allocation is long-only. The ordinary case is projected onto one
  unit simplex; taxable and tax-advantaged allocations are projected onto two
  simplices with their prescribed account totals.
- The total-net-worth optimizer retains the upstream treatment of human
  capital and liabilities as asset-exposure vectors.
- Assets with zero volatility remain deterministic in simulated returns.

## Differences from R

The upstream optimizer delegates to `nloptr`. The Fortran implementation uses
analytic first derivatives of expected CRRA utility and projected backtracking.
Its benchmark results match the upstream tests, but iteration paths and status
codes are necessarily different.

The R package evaluates arbitrary tidyverse formulas against nested household
records. General R expression evaluation is not appropriate in a compiled
Fortran library. Procedure callbacks and explicit cash-flow vectors provide the
same computational flexibility while remaining type safe.

The upstream `simulate_single_scenario()` uses the initial taxable and
 tax-advantaged weights when computing the certainty-equivalent discount rate
for discretionary spending. The Fortran simulator preserves that convention.

The source uses internal procedure callbacks for native optimizers. Some GNU
linkers on Unix may print an executable-stack notice for such callbacks; this
is a compiler/linker implementation detail and does not affect numerical
results. Windows GNU builds generally do not emit that notice.

## Data

The original binary `life_tables.rda` is retained only in the upstream
snapshot. It is not decoded or redistributed as generated Fortran source.
Applications should read their preferred mortality data and supply vectors of
ages and one-year mortality rates to `fit_gompertz_mortality()`.
