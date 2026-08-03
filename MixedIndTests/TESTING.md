# Testing

With FPM:

```text
fpm test
fpm run
fpm run --example independence_example
```

With GNU Fortran:

```text
./scripts/test_gfortran.sh
```

The strict script enables Fortran 2018 conformance, all common warnings,
warnings as errors, bounds and allocation checking, and traps for invalid,
divide-by-zero, and overflow exceptions.

The four test programs cover:

1. exact original-C regression values for data preparation, pair dependence,
   nonserial statistics, multiplier matrices, and Moebius scores;
2. exact scalar/vector serial regression values and bootstrap calculations;
3. high-level pairwise, serial, multivariate, bootstrap, and Moebius tests;
4. all simulation families, all margin categories, reproducibility, and
   data-driven lag selection.
