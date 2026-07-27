# Validation

## Recorded environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- GNU Make
- Runtime checking and backtraces enabled

Checked flags:

```text
-std=f2018 -O0 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace
```

`fpm.toml` is supplied, but `fpm` was not installed in the validation
environment and therefore is not claimed as tested.

## Reproduce the checked build

```text
make clean
make check
```

`make check` first verifies GPL-3.0-only headers in every Fortran source file,
then compiles all modules and runs these programs.

### `test_dcc`

- Gaussian DCC simulation, filtering, constrained fitting, and forecasting
- Student ADCC simulation and filtering
- Correlation diagonals and bounds
- Nonnegative coefficients and stationarity

### `test_core`

- Gaussian copula identity density
- FDCC filtering
- GO-GARCH covariance and volatility construction
- VAR fitting, filtering, forecasting, and simulation
- FastICA execution
- DCC constancy diagnostic
- Mardia skewness and kurtosis diagnostics

### `test_extended`

- Multivariate Normal, standardized Student-t, and Laplace densities/simulation
- Student and Laplace scalar CDF/quantile transforms
- Weighted margins and weighted-margin paths
- General DCC(2,1), history-aware forecasting, Laplace DCC fit, and Student
  shape estimation
- Grouped FDCC fitting, forecasting, and simulation
- Static/dynamic Gaussian and Student copula fitting, likelihoods,
  transformations, and simulation
- FastICA and RADICAL reconstruction
- GO-GARCH moments, betas, fitted model, filtering, forecasting, simulation,
  and rolling forecasts
- Raw two-step DCC and raw-return rolling covariance forecasts
- Model-level dynamic Gaussian Copula-GARCH fit/filter/simulation
- FFT Normal-grid convolution with mean, variance, density, CDF, and quantile
  checks
- Conditional DCC scenarios, higher moments, and portfolio projection
- Robust VARX with an exogenous regressor

Expected output:

```text
DCC tests passed.
Core tests passed.
Extended computational tests passed.
```

## Additional release checks

The release procedure also performs:

```text
make clean
make all fit_csv
./build/demo_rmgarch
./build/fit_csv data/sample_returns.csv
```

and an optimized warning-enabled rebuild with:

```text
make clean
make FFLAGS="-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface" all fit_csv
```

These checks establish internal execution and invariant consistency. They do
not establish exhaustive numerical equivalence with R `rmgarch` 1.4-2.
