# Testing

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- GNU/Linux x86-64
- Python 3 with NumPy and SciPy for independent reference values
- No `fpm` executable was installed in the translation environment

## Strict debug build

The library and tests were compiled with:

```text
-std=f2018
-Wall -Wextra -Werror -pedantic
-fcheck=all
-ffpe-trap=invalid,zero,overflow
-O0 -g
```

## Optimized build

A separate clean build used:

```text
-std=f2018
-Wall -Wextra -Werror -pedantic
-O3
```

## Test coverage

`test/test_risk.f90` checks:

- standard-normal density, CDF, and quantile references;
- normal VaR and ES analytical values;
- agreement of quantile-, CDF-, and PDF-based ES;
- location-scale risk;
- log-return to simple-return VaR and ES;
- empirical R type-7 VaR;
- empirical ES tail selection;
- scalar, vector, and matrix overloads;
- arbitrary logistic callback functions;
- Student-t and standardized-t SciPy references;
- GED SciPy references;
- invalid-probability status handling.

`test/test_garch.f90` checks:

- model validation and unconditional variance;
- exact observation and variance recursions;
- deterministic repeatability for a fixed seed;
- normal innovation sample moments;
- one-step and multi-step variance forecasts;
- plug-in predictive intervals;
- Monte Carlo interval ordering;
- standardized-t innovation variance;
- GED innovation variance;
- invalid model rejection.

## Independent references

Fixed references were generated independently with SciPy 1.17.0 for:

- normal quantiles and expected shortfall;
- Student-t and standardized-t quantiles;
- generalized-error quantiles and density;
- transformed lognormal lower-tail expected shortfall.

The original R test files are retained under `original/tests/testthat/`. R was not installed in the build environment, so the `.RDS` random-stream fixtures were not executed. Exact stream equality would not be expected because this port deliberately uses a portable explicit RNG rather than R's RNG.

## FPM manifest

`fpm.toml` was parsed with Python's TOML parser and all application, example, and test targets were compiled directly with `gfortran`. On a system with FPM, use:

```text
fpm test
```
