# REN modern Fortran

A modern Fortran/FPM translation of the computational code in REN 0.1.0,
"Regularization Ensemble for Robust Portfolio Optimization." Plotting code and
R-only parallel/runtime infrastructure are omitted.

## Build

```text
fpm build
fpm test
fpm run --target demo_ren
```

The package is self-contained. The compatible supplied corpcor translation is
vendored as a local FPM dependency. REN's required Gaussian LASSO, ridge, and
elastic-net routines are implemented directly in this package.

## Main API

```fortran
use ren, only : dp, analysis_options, analysis_result, perform_analysis

type(analysis_options) :: options
type(analysis_result) :: result

options%cluster_repetitions = 100
options%stochastic_samples = 1000
call perform_analysis(returns_percent, month_index, yyyymmdd, result, options)
```

`returns_percent(row, asset)` follows the R package convention: a value of
`1.25` means a 1.25% return. The output contains monthly weights and turnover,
daily gross returns, arithmetic cumulative returns, compounded wealth indices,
Sharpe ratios, annualized volatility, and standard maximum drawdown. No plotting
objects are produced.

The individual translated constructors are also public: `po_cols`, `po_jm`,
`po_avg`, `po_gross_exp`, `po_cov_shrink`, `buh_clust`, `po_bhu`, `po_tzt`,
`po_sw`, and `po_sw_lasso`.

See `API_MAP.md` and `PORTING_NOTES.md` for exact scope and source differences.
