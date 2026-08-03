# Testing

## FPM

```text
fpm test
```

## GNU Fortran strict build

```text
./scripts/test_gfortran.sh
```

This uses Fortran 2018 conformance, warnings as errors, array/runtime checking,
floating-point traps, and backtraces.

## GNU Fortran optimized build

```text
./scripts/test_gfortran_optimized.sh
```

This uses `-O3` and warnings as errors.

The four test programs cover:

1. Scalar and column-wise matrix moments, skewness, and kurtosis.
2. Raw/central conversions and corrected/legacy cumulants.
3. Regression values for all four hypothesis tests and alternatives.
4. NaN, constant-sample, insufficient-sample, and degenerate-data paths.
