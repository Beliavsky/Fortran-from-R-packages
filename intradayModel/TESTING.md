# Testing

Four independent test programs are included:

1. `test_kalman`: filtering, smoothing, covariance symmetry, and state-error reduction.
2. `test_fit`: ordinary EM convergence and recovery of synthetic model parameters.
3. `test_fixed_acceleration`: accelerated EM, fixed parameters, initialization, and the
   one-iteration nonconvergence path.
4. `test_decompose`: analysis/forecast components, burn-in, wrapper behavior, error
   metrics, and invalid-day removal.

The Unix validation scripts compile with GNU Fortran 14 using either:

```text
-std=f2018 -Wall -Wextra -Wpedantic -Werror
-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
```

or an optimized configuration using `-O3 -Werror`.
