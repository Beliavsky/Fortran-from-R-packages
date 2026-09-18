# Validation

## Direct GNU Fortran validation

The maintained `varmapack` Fortran sources, deterministic tests, and example
were compiled with GNU Fortran 14.2.0 using:

```text
-std=f2018 -O2 -fcheck=all -Wall -Wextra -Werror=line-truncation
```

The actual translated `randompack` source was used. Because the target
repository's `rfortran-linalg` checkout was not available inside this execution
container, direct validation used an **external validation-only shim** exposing
the current `r_linalg` high-level interfaces required by `varmapack` and
`randompack`. That shim used the container's LAPACK/BLAS libraries solely for
validation. It is outside the package tree and is not present in the ZIP.
The package manifest itself contains no system BLAS/LAPACK links and depends on
`../rfortran-linalg` as required.

The deterministic test executable completed with:

```text
All varmapack tests passed.
```

The example completed with output of the form:

```text
spectral radius:   0.5000
theoretical covariances:    1.33333   0.66667   0.33333
first simulated path:    1.08988  -0.41079   1.27199  -1.62635
```

The suite covers, among other paths:

- AR(1) spectral radius, ACVF, PSI, and IRF;
- multivariate PSI matrix orientation;
- ML and corrected sample autocovariance normalization;
- covariance-to-correlation conversion;
- all testcase metadata plus named/indexed, random-repeatability, and
  rho-controlled construction;
- white noise and singular innovation covariance;
- stationary seeded VARMA repeatability;
- time-dependent mean paths;
- nonstationary supplied-start VARMA recursion;
- VARMAX startup conditioning and exogenous recurrence;
- zero-start VARMAX simulation.

## FPM validation limitation

`fpm` is not installed in this execution environment. The required commands
were explicitly attempted before packaging:

```text
fpm build
fpm test
fpm run --example basic_varmapack
fpm clean --all
```

and could not execute because the `fpm` command is unavailable. Therefore this
document does **not** claim a successful FPM integration build. A final
release-quality check in the target repository should place `varmapack`,
`randompack`, and `rfortran-linalg` as sibling directories and run those four
commands on Windows/gfortran as requested.

## Source/package audit

Before the archive is created, automated checks are run for:

- maintained Fortran line lengths;
- free-form `.f90` source only;
- one declaration per dummy argument with explicit `INTENT`/`VALUE` and a
  same-line trailing FORD `!!` description;
- forbidden `double precision`, `real*8`, `kind(0.0d0)`, and `d0`-style
  literals;
- executable semicolon-separated statements;
- self-comparison NaN idioms;
- duplicate maintained Fortran source files;
- coverage metadata/path consistency;
- absence of object/module/executable/cache/ZIP build products from the
  package tree;
- absence of vendored sibling dependency source.
