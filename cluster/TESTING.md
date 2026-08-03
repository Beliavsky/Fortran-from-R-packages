# Testing

Four test programs cover:

1. Numeric and mixed `daisy`, PAM, CLARA, and FANNY.
2. AGNES, DIANA, MONA, hierarchy dimensions, ordering, and coefficients.
3. Silhouette widths, medoids, distance-size and triangular-index helpers,
   `max_se`, and gap statistics through a typed callback.
4. Minimum-volume enclosing ellipsoid construction, containment, volume, and
   generated boundary points.

Validation commands:

```text
scripts/test_gfortran.sh
scripts/test_gfortran_optimized.sh
```

The strict configuration enables Fortran 2018 conformance, warnings as errors,
bounds checking, and floating-point traps. The optimized configuration uses
`-O3` with warnings as errors.
