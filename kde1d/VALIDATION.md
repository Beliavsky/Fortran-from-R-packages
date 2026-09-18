# Validation

Final validation date: 2026-09-15.

## FPM commands

The requested FPM commands were attempted from the package root:

```text
fpm build
fpm test
fpm run --example basic_kde1d
fpm clean --all
```

They could not be executed in this sandbox because no `fpm` executable is
installed; each command returned status 127 (`fpm: command not found`). This is
an environment limitation, not a claimed successful FPM run. The manifest was
parsed separately and the sibling `rfortran-core` dependency/API names were
checked against the current `Fortran-from-R-packages` repository.

## Direct GNU Fortran validation

GNU Fortran 14.2.0 successfully compiled and ran the maintained package source,
test program, and example using API-compatible validation modules for the
external `rfortran-core` procedures. The package source was compiled with:

```text
-std=f2018 -Wall -Wextra -Werror=line-truncation -fcheck=all -fbacktrace -O0
```

The deterministic test executable reported:

```text
All kde1d tests passed.
```

The example also ran successfully. The test suite exercises fixed and automatic
bandwidths, local-polynomial degrees 0/1/2, unbounded and one-/two-sided finite
supports, endpoint repair, scale equivariance of the one-sided transform,
weighted fits, discrete normalization/CDF/quantiles, zero-inflated fits
(including bounded and all-zero cases), density/CDF/quantile inversion, and
pseudo/quasi sampling paths.

Compiler warnings from exact real comparisons are expected in a few places
where equality is part of the upstream semantics (integer-level tests, zero
atoms, equal weights, and ties). No fast-math or finite-math-only options were
used.

## Static release audit

The final maintained tree passed checks for:

- free-form `.f90` translated/new source only;
- a single `real(dp)` kind imported from `rfortran-core`;
- no `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals;
- no semicolon-separated Fortran statements;
- no self-comparison NaN tests;
- explicit `INTENT` or `VALUE` for every dummy argument;
- one declaration per dummy argument with a meaningful trailing FORD `!!`
  comment;
- no maintained Fortran line longer than 132 characters;
- no duplicate maintained Fortran source files;
- no copied `rfortran-core`, BLAS, LAPACK, ARPACK, `r.f90`, or `r_mod.f90`
  source;
- no system BLAS/LAPACK links or `rfortran-compat` dependency;
- all eight coverage mappings and source paths present and internally
  consistent;
- no object files, module files, executables, caches, libraries, or ZIP files
  inside the package tree.

Because FPM is unavailable, `fpm clean --all` could not perform cleanup. The
package tree was instead inspected directly for build products before archive
creation.
