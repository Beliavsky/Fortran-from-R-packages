# Validation

Validation was performed with GNU Fortran 14.2.0 and system BLAS/LAPACK.

Compiler checks used for the root translation:

```text
-std=f2018 -fcheck=all -Werror=implicit-interface
```

No unlimited free-form line-length option is required; root source, tests, and example
lines are kept within standard free-form limits.

The current regression suite contains nine programs:

1. `test_core` — transforms/features/benchmarks/Croston/DM checks.
2. `test_ets_arima` — ETS, ARIMA, forecasting, and automatic order selection.
3. `test_dependencies` — `urca`-backed differencing, `fracdiff` ARFIMA, and `nnet` nnetar.
4. `test_bats` — BATS/TBATS state placement and high-level fitting/forecasting.
5. `test_helpers` — DSHW, decomposition, cleaning, and bootstrap.
6. `test_regression_spline` — regression/modelAR/spline fitting and forecasting.
7. `test_cv_bagging` — rolling-origin CV and bagging.
8. `test_parity_v02` — stationary ARMA Gaussian likelihood, OCSB critical values,
   fixed-parameter ETS and variance recursion, generic bagging/CVar, STL seasonal strength,
   STLF reseasonalization, TBATS harmonic search, calendar helpers, multi-season Fourier
   handling, and tapered ACF.
9. `test_parity_v03` — ARIMAX and structural refits, `auto_arima` truncation, closed-form
   diffuse random-walk likelihood, missing-observation ML, integrated simulation, ETS
   class-3 and simulated intervals, type-8 quantiles, tapered bootstrap CIs, and BATS/TBATS
   forced-control/refit behavior.

All nine programs pass under the flags above. The example program is also compiled and executed in the development validation workflow.
FPM itself was not available in that environment; `fpm.toml` is retained for normal FPM
builds and the equivalent package-local dependency graph was compiled manually.
