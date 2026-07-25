# Validation

The release is validated with GNU Fortran 14.2.0, LAPACK, and BLAS.

## Runtime-checked build

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

Executed by `make check`:

- License audit over every Fortran source file
- `test_core`
- `test_fit`
- `test_inference`
- `demo_garchx`
- `sparse_garchx_example`
- `fit_csv` using ordinary covariance
- `fit_csv` using robust covariance
- `fit_csv` using HAC covariance

## Optimized build

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
```

Executed by `make optimized-check`, which rebuilds and reruns the same checks.

## Numerical checks

`test_core` checks:

- Sparse selected-lag recursion against fixed independent values
- Exact recursive variance derivatives against central finite differences
- Objective value against a direct calculation
- Zero-return objective mode
- Simulation identity with supplied innovations
- Vector and matrix lag/difference operations
- Normal and Student-t probability/quantile inversions

`test_fit` checks:

- Simulation and fitting of asymmetric GARCH-X with a covariate
- Positive fitted variances and finite objective/Hessian output
- Ordinary, robust, and HAC covariance symmetry and finiteness
- Multi-step bootstrap forecasts and path averaging
- Empirical conditional-quantile paths
- Fixed-parameter and re-estimated refits
- Simulation-based ordinary asymptotic covariance

`test_inference` checks:

- Multivariate Normal sample means and covariance
- Student-t confidence intervals
- Boundary-null t statistic and one-sided p-value
- Simulated boundary-Wald critical values and their ordering

No exact R optimizer or RNG equivalence is claimed.
