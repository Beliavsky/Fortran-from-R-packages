# Testing

## FPM

```console
fpm build
fpm test
fpm run
fpm run --example time_series_diagnostics
fpm run --example arima_example
fpm run --example apca_example
```

## GNU Fortran scripts

Linux and macOS with GNU Fortran:

```console
./run_gfortran_tests.sh
```

Windows command prompt with GNU Fortran on `PATH`:

```console
run_gfortran_tests.bat
```

The shell script defaults to strict runtime-checked flags:

```text
-std=f2018 -Wall -Wextra -Werror -pedantic
-fcheck=all -ffpe-trap=invalid,zero,overflow
-ffree-line-length-none -O0 -g
```

An optimized validation can be run with:

```console
FFLAGS='-std=f2018 -Wall -Wextra -Werror -pedantic -ffree-line-length-none -O3' \
  ./run_gfortran_tests.sh
```

## Test programs

- `test_finance_stats`: moments, interest, return conversion, month conversion, and conjugate roots
- `test_time_series`: ACF, covariance, PACF, cross-ACF, Box tests, ARCH LM, ARMA ACF, roots, and stationarity
- `test_apca`: one-factor APCA dimensions, eigenvalue ordering, normalization, and explanatory power
- `test_arima`: AR fitting, regressors, differencing, seasonal AR fitting, variance, and residual diagnostics
