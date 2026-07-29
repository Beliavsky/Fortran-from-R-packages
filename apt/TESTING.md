# Testing

Four independent test programs are included.

- `test_statistics`: distribution functions and OLS fixed references.
- `test_reference`: fixed NumPy/SciPy references for TAR, MTAR, symmetric ECM, and asymmetric ECM coefficients, SSE, AIC, and F statistics.
- `test_search`: common-window lag selection and trimmed threshold-search invariants.
- `test_ecm_tests`: H1/H2 restriction references, Durbin-Watson statistics, and Ljung-Box probabilities.

The reference series are generated deterministically in both Python and Fortran. NumPy was used for matrix calculations and SciPy for Student-t, F, and chi-square probabilities. The generation method is recorded in `REFERENCE_GENERATION.md`.

The release script builds two configurations:

1. strict debug: Fortran 2018, warnings as errors, bounds/runtime checks, and floating-point traps;
2. optimized: `-O3` with warnings as errors.

Every test, application, and example is run in both configurations.
