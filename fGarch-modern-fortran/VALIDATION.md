# Validation status

## Completed checks

The project was compiled with GNU Fortran 14.2 using Fortran 2018 mode, warnings, explicit-interface diagnostics, bounds/runtime checks, and backtraces.

The included tests check:

- Normal density at zero.
- GED with shape 2 reducing to the standard normal density.
- Standardized Student-t density at zero.
- Student-t and GED CDF/quantile round trips.
- Monotonic corrected skew-normal, skew-t, and skew-GED quantiles around probabilities 0.49, 0.50, and 0.51.
- Skew distribution CDF/quantile round trips.
- Unit second absolute moments for standardized normal, Student-t, and GED.
- GARCH simulation, filtering, likelihood, constrained fitting, and forecasts.
- APARCH simulation and likelihood under skew GED innovations.
- Normal GARCH kappa and persistence identities.
- The no-argument public demonstration.
- The numeric CSV fitting example.

## Commands used

```text
make test
make all
./build/demo_fgarch
./build/fit_csv data/sample_returns.csv
```

## Not yet validated comprehensively

- Exact equivalence to all R fGarch optimizers and initialization paths.
- All package unit tests and historical bug cases.
- Hessian, covariance matrix, and standard-error calculations.
- ARMA parameter estimation.
- General p/q parameter estimation.
- EGARCH estimation and multi-step EGARCH forecasting.
- Very heavy-tailed edge cases near nonexistent moments.
- Cross-compiler and cross-platform behavior.
- Reproduction of R random-number streams.

This is an experimental release, not a claim of full numerical equivalence.
