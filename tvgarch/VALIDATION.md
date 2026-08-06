# Validation

The test suite covers:

1. Logistic transition identities, component construction, parameter packing,
   and binary combination enumeration.
2. Univariate TV-GARCH simulation identities and positive variances.
3. Simulated TV-GARCH fitting, location/size recovery, forecasts, and quantile paths.
4. CCC/DCC multivariate simulation, DCC estimation, marginal fitting, lower-
   triangle extraction, and cross-variance spillover fitting.
5. Robust and nonrobust transition-order test statistics and valid p-values.

Commands executed in the validation environment:

```text
make check
make optimized-check
```

Compiler: GNU Fortran 14.2.0. BLAS/LAPACK were linked from the system.
All test executables and the demonstration completed successfully in checked
and optimized configurations.
