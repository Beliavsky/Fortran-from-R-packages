# Validation

Validation was performed on the translated source with GNU Fortran 14.2.0 using strict Fortran 2018 mode, warnings, line-truncation errors, optimization, and runtime checking.

Because this execution environment does not contain the Fortran Package Manager and does not contain the sibling repository checkouts as an FPM workspace, the source was compiled against small validation-only interface shims matching the public APIs used from `rfortran-core::r_kinds` and `rfortran-linalg::solve_system`. Those shims are not included in this package ZIP.

The deterministic test suite exercises all 18 mapped APIs, including weighted windows, missing logical values, type-2 quantiles, unbiased moments, complete-row matrix handling, pairwise cubes, and rolling linear regression. It reports:

```text
All roll tests passed.
```

The example also compiles and runs successfully.

The following required release commands were explicitly attempted in the package directory:

```text
fpm build
fpm test
fpm run --example basic_roll
fpm clean --all
```

They could not run because `fpm` is not installed in this environment. Therefore this document does **not** claim a completed FPM integration build against the real sibling dependencies. A release checkout should rerun those commands with `roll`, `rfortran-core`, and `rfortran-linalg` as sibling top-level directories.

Static checks performed before packaging verify free-form source only, the shared `dp` kind, no disallowed semicolon-separated statements, no self-comparison NaN idioms, no forbidden real-kind forms, no system BLAS/LAPACK links, no copied dependency source, no duplicate maintained Fortran files, trailing FORD comments on dummy-argument declarations, and absence of build products/caches in the archive.
