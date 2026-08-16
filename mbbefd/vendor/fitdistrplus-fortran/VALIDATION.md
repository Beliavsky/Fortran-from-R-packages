# Validation

Five deterministic test programs cover:

1. Distribution densities, CDFs, inverse CDFs, moments, and factory aliases.
2. MLE, MME, QME, MGE, MSE, automatic starts, and MLE covariance.
3. Mixed exact, left-, right-, and interval-censored exponential fitting.
4. Descriptive statistics, KS/CvM/AD diagnostics, grouped chi-square, and parameter bounds.
5. Parametric bootstrap, CDF bands, custom callback fitting, and interval-censoring NPMLE.

The checked GNU Fortran build uses:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -pedantic -fcheck=all -fbacktrace
```

The optimized build uses `-O3` and the same strict diagnostics, except for
GNU's optimizer-only `uninitialized` and `maybe-uninitialized` analysis of
allocatable descriptors in local derived-type results.

The test executables have a non-executable GNU stack (`GNU_STACK RW`).
