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

The tests cover:

1. Deterministic zero-volatility diffusion values in corrected and legacy mode.
2. Fixed-seed reproducibility, stored paths, discounting, and standard errors.
3. Equality of `jdm_bs(lambda=0)` and `normal_bs` under the same random stream.
4. Pathwise Monte Carlo put-call identities.
5. Shared-event propagation and expected Poisson event counts.
6. Zero transmission, invalid dimensions, negative parameters, and NaN inputs.

The strict script uses Fortran 2018, all common warnings as errors, bounds and
runtime checks, floating-point traps, and backtraces. The optimized script uses
`-O3` with warnings as errors.
