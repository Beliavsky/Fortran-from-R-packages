# Testing

## FPM

```text
fpm test
```

## GNU Fortran

Unix-like systems:

```text
./scripts/test_gfortran.sh
./scripts/test_gfortran_optimized.sh
```

Windows with GNU Fortran tools on `PATH`:

```text
scripts\test_gfortran.bat
```

## Test coverage

- Fixed BNS, Amin, and Amed regression values.
- Scalar and period jump-test p-values.
- Benjamini-Hochberg adjustment and p-value combination.
- Independent and dependent Fisher/Stouffer pooling.
- Minimum and maximum p-value pooling.
- Exact caller-supplied-shock recursions for all three native model kernels.
- Seed reproducibility and output dimensions for all five simulators.
- No-jump limits and nonnegative CIR variance paths.
- Invalid methods, dimensions, correlations, covariance matrices, and
  degenerate return samples.

The strict script enables bounds checking, floating-point traps, backtraces,
and warnings as errors. The optimized script uses `-O3` with warnings as
errors.
