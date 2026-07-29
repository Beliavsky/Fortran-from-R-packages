# Testing

## One-command test

```bash
./scripts/test_gfortran.sh all
```

This creates clean `debug` and `release` build directories.

Debug flags include:

```text
-std=f2018 -O0 -g -Wall -Wextra -Werror
-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
```

Release flags include:

```text
-std=f2018 -O3 -Wall -Wextra -Werror
```

BLAS and LAPACK are linked explicitly.

## Test programs

### `test_core`

- fixed independent NumPy likelihood references;
- full, diagonal, scalar, and asymmetric branches;
- model-validity checks;
- fixed-innovation simulation references.

Reference log likelihoods:

```text
full                  7.559158365670399
asymmetric full       7.555407727379082
diagonal              7.559188653634845
scalar                5.560859831717117
asymmetric scalar     6.051637028692405
```

The reference generator is `scripts/generate_references.py`.

### `test_matrix`

- lower-triangle vectorization;
- elimination and duplication identities;
- commutation identity;
- lag-matrix construction.

### `test_inference`

- deterministic and randomized starting values;
- per-observation score versus independent finite differences;
- Hessian symmetry;
- OPG and sandwich covariance symmetry;
- public fixed-innovation simulation;
- Moore-Penrose generalized-inverse recovery.

### `test_workflow`

- asymmetric scalar simulation and fitting;
- stationarity and fitted covariance paths;
- multi-step forecasting;
- normal portfolio VaR;
- Kupiec/Christoffersen backtesting;
- VIRF generation and delta-method intervals;
- portmanteau diagnostics;
- rolling-window backtesting;
- Monte Carlo parameter-recovery evaluation.

## Runnable targets

The test script also builds and executes:

- `app/bekks_demo.f90`
- `example/asymmetric_bekk.f90`
- `example/portfolio_var.f90`

## FPM

The repository follows standard FPM layout. With FPM installed:

```bash
fpm test
fpm run
fpm run --example asymmetric_bekk
fpm run --example portfolio_var
```
