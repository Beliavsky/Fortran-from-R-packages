# Testing

Run with FPM:

```text
fpm test
```

Or with GNU Fortran:

```text
./scripts/test_gfortran.sh
./scripts/test_gfortran_optimized.sh
```

The four tests cover:

- Normal and Laplace special cases
- Reference values for shapes 0.5, 1.3, 2, and 5
- Density/CDF logarithmic and tail options
- Quantile/CDF round trips down to probability `1e-8`
- Upper-tail log-probability semantics
- Gamma-function identities and inverse calculations
- Seeded random reproducibility and sample moments
- Invalid parameters and endpoint probabilities

The strict script enables bounds checking, floating-point traps, backtraces,
and warnings as errors. The optimized script uses `-O3 -Werror`.
