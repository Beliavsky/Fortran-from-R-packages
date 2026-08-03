# esback-fortran

Modern Fortran 2018/FPM translation of the computational algorithms in R package **esback 0.3.1**.

The library implements:

- McNeil-Frey exceedance-residual backtests, with simple and volatility-standardized residuals
- Nolde-Ziegel conditional-calibration backtests, including simple/general instruments and Hommel or Bonferroni correction
- Bayer-Dimitriadis expected-shortfall-regression backtests:
  - strict ESR
  - auxiliary ESR
  - strict-intercept ESR
- asymptotic and iid-bootstrap p-values
- the required `esreg` numerical subset: joint VaR/ES Fissler-Ziegel loss, quantile-regression starts, iterated local search, misspecification-robust sandwich covariance, density-at-quantile estimation, conditional location-scale fitting, and truncated conditional variance

## Build

```sh
fpm test
fpm run
```

The library links LAPACK and BLAS. A reproducible GNU Fortran build is also provided:

```sh
./tools/test_gfortran.sh strict
./tools/test_gfortran.sh optimized
```

## Minimal use

```fortran
use esback

type(er_backtest_result) :: er
call er_backtest(returns, var_forecasts, es_forecasts, er, &
  volatility_forecasts, 1000)
```

See `example/`, `API.md`, and `PORTING.md`.

## Scope

R data frames, formulas, S3 methods, printing, and the bundled serialized dataset are not runtime dependencies. Numerical arrays and typed result structures replace them. The original source tree is retained under `original/` for provenance.

## License

GPL-3.0-only, matching esback and the embedded esreg-derived numerical subset. See `LICENSE` and `NOTICE.md`.
