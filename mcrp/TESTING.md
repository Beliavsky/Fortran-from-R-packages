# Testing

## Build configurations

The release was compiled and run in two configurations with GNU Fortran 14.2:

### Strict

```sh
./build_gfortran.sh strict
```

Uses Fortran 2018, warnings as errors, runtime bounds and consistency checks,
backtraces, and floating-point traps for invalid operations, division by zero,
and overflow.

### Optimized

```sh
./build_gfortran.sh release
```

Uses `-O3` with the same language and warning checks.

## Test programs

### `test_moments`

Checks:

- dimensions and independent reference entries for `M2`, `M3`, and `M4`
- exact `pm2`, `pm3`, and `pm4` contractions
- derivative reference values
- raw and standardized portfolio measures
- contribution identities for percentage and absolute modes
- all original alias families

Reference values were generated independently with NumPy using explicit tensor
loops. The script is `tools/generate_reference.py`.

### `test_optimizer`

Checks:

- objective scale invariance
- nonnegative bounded optimization
- large reduction in the pure variance-risk-parity objective
- nearly equal variance contributions
- reduction of the combined variance/skewness/kurtosis objective
- contribution-sum identities at the fitted portfolio
- L1 normalization of final weights

## Runnable targets

Both builds also compile and execute:

- `mcrp_demo`
- `moment_example`

## Results

All strict and optimized tests and runnable targets passed. The final release
archive was extracted and rebuilt from its contents before publication.
