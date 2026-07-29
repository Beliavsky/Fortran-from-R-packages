# Testing

The package was tested with GNU Fortran 14.2.0 in two configurations.

## Strict debug configuration

```console
gfortran -std=f2018 -Wall -Wextra -Werror -pedantic \
  -fcheck=all -ffpe-trap=invalid,zero,overflow \
  -ffree-line-length-none -O0 -g
```

## Optimized configuration

```console
gfortran -std=f2018 -Wall -Wextra -Werror -pedantic \
  -ffree-line-length-none -O3
```

## Test programs

### `test_smoothing`

- independent fixed references for left-boundary, interior, and right-boundary local-polynomial estimates
- independent equivalent-kernel weight reference
- derivative scaling
- exact cubic-polynomial reproduction
- independent Nadaraya-Watson references

### `test_cf0_arma`

- independent Bühlmann lag-window value and selected lags
- seeded ARMA(1,1) simulation and parameter recovery
- MA-infinity recursion
- complete information-criterion matrix construction and order selection

### `test_estimation`

- all ten named `msmooth` algorithms
- first- and second-derivative iterative plug-in estimators
- confidence-bound ordering
- trend/residual reconstruction identity
- valid automatic bandwidth ranges

### `test_forecast`

- normal ARMA forecast intervals
- deterministic forward-bootstrap intervals and exported error matrix
- combined trend/residual forecasting
- rolling holdout forecasts, interval breaches, MASE, and RMSSE

All applications and examples are compiled and executed as part of the clean-build script.

## Running

With FPM:

```console
fpm test
```

Without FPM:

```console
./run_gfortran_tests.sh
```
