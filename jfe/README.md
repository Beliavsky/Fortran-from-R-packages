# JFE-fortran

Modern Fortran translation of the computational core of the R package
**JFE 2.5.11**.

The library provides performance and risk measures for numeric return vectors,
CAPM-relative measures, drawdown statistics, tail-risk ratios, annualized
summaries, and the Durbin h test. It is a self-contained Fortran Package
Manager project and has no external numerical dependencies.

## Build with FPM

```text
fpm build
fpm test
fpm run demo_jfe
fpm run --example performance_indices
```

GNU Fortran validation scripts are also supplied:

```text
./run_tests.sh
./run_release_tests.sh
```

## Interface conventions

R names containing dots are mapped to underscores. For example:

- `Return.annualized` -> `return_annualized`
- `SharpeRatio.annualized` -> `sharpe_ratio_annualized`
- `CAPM.jensenAlpha` -> `capm_jensen_alpha`
- `table.AnnualizedReturns` -> `table_annualized_returns`

R `xts` columns are represented by one-dimensional `real(dp)` arrays. Apply a
procedure column by column for a matrix of return series. IEEE NaNs are omitted
where the R implementation used `na.omit` or `na.rm=TRUE`.

Frequencies are explicit through a `scale` argument: commonly 252 for daily,
52 for weekly, 12 for monthly, 4 for quarterly, and 1 for annual data.

## Scope

All 32 non-download exported computational routines are represented. Network
download functions (`getEER`, `getFed`, `getFrench.Factors`, and
`getFrench.Portfolios`) and R-specific `xts`, plotting, data-frame, and model
object infrastructure are not part of the compiled library.

See `API.md`, `PORTING.md`, and `TRANSLATION_COVERAGE.md` for details.

## License

GPL-2.0-or-later, matching the upstream package. Original package sources and
metadata are retained under `original/JFE/`.
