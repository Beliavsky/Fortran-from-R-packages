# Testing

Run all FPM tests with:

```text
fpm test
```

For stricter direct GNU Fortran validation:

```text
sh scripts/test_gfortran.sh
sh scripts/test_gfortran_optimized.sh
```

On Windows with GNU Fortran available in `PATH`:

```text
scripts\test_gfortran.bat
```

The permanent tests cover:

- weighted OLS quantities, generic meat, HC0/HC3, and OPG covariance;
- all representative kernel formulas, PAVA, HAC lag algebra, prewhitening,
  and automatic bandwidth selection;
- one-way and two-way cluster inclusion-exclusion;
- aggregate and within-cluster panel longitudinal covariance;
- panel-corrected covariance;
- bootstrap and jackknife scaling, deterministic wild bootstrap, symmetry,
  and PSD repair.

Small test matrices use independently calculated expected values rather than
self-comparison with another path through the library.
