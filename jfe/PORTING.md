# Porting notes

## R objects and time indexes

The R package accepts `xts` objects, infers periodicity from their indexes, and
returns matrices with row and column labels. The Fortran API accepts plain
numeric vectors and requires annualization scale explicitly. This avoids hidden
calendar assumptions and keeps the library independent of a date-time package.

## Missing values

IEEE NaNs are removed in routines corresponding to upstream `na.omit` or
`na.rm=TRUE` behavior. Paired calculations remove an observation whenever
either series is nonfinite.

## Return-unit conventions

The upstream source uses two different conventions:

- `Return.annualized`, standard drawdowns, and most ratios use decimal returns.
- `DrawdownPeak`, `UlcerIndex`, `PainIndex`, and Burke drawdown episodes divide
  returns by 100 and therefore treat inputs as percentage returns.

The translation preserves this distinction for numerical compatibility rather
than silently normalizing units.

## Corrections and safeguards

- Upstream `BurkeRatio(..., modified=TRUE)` reads `result` before assigning it.
  The Fortran routine implements the standard modified Burke denominator,
  equivalent to multiplying the ordinary ratio by `sqrt(n)`.
- The nonannualized upstream `SharpeRatio` accepts `Rf` but does not subtract it.
  The Fortran procedure subtracts `rf`, matching the documented definition.
- Invalid probabilities, insufficient observations, zero denominators, and an
  invalid Durbin-h square-root term return IEEE NaN plus an explicit status.
- Historical VaR follows R quantile `type=1`, as in the original package.

## Omitted functions

`getEER`, `getFed`, `getFrench.Factors`, and `getFrench.Portfolios` perform
network access and file parsing. They are retained in `original/JFE/R/` but are
not compiled. The translation also omits plotting, labels, formula/model-object
introspection, and bundled `.rda` data loading.
