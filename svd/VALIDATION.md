# Validation

This record distinguishes checks actually executed in the translation environment from
FPM/Windows checks that require the target repository layout.

## Translation environment

- GNU Fortran 14.2.0, Linux x86-64.
- `fpm` was not installed in the environment.
- `fprettify` was not installed in the environment.
- Network access from the compiler container was unavailable, so sibling dependencies
  could not be cloned into the compiler filesystem.

The current public interfaces of sibling `RSpectra` and `rfortran-linalg` were inspected
before implementation. The maintained `svd` source was then compiled against small
API-compatible validation stand-ins kept **outside** the release tree. Those stand-ins use
the host LAPACK only to exercise the adapter layer during this local validation; they are
not package dependencies, are not included in the ZIP, and do not change `fpm.toml`.
The released package itself has no system BLAS/LAPACK/ARPACK link directive.

This adapter validation is strong evidence for the maintained Fortran source and its
numerical mappings, but it is not represented as a successful full FPM integration build
against the real sibling source trees.

## Checked GNU Fortran build

All maintained modules, tests, the demo, and examples were compiled from a clean build
with:

```text
-std=f2018 -O0 -g -Wall -Wextra -Werror -Wconversion-extra
-Wimplicit-interface -fcheck=all -fbacktrace -ffree-line-length-none
```

Results:

```text
test_extmat: PASS
test_propack: PASS
test_trlan: PASS
test_complex: PASS
svd_demo: PASS
dense_svd: PASS
matrix_free_svd: PASS
```

## Optimized GNU Fortran build

The complete maintained source graph was rebuilt from scratch with:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wconversion-extra
-Wimplicit-interface -ffree-line-length-none
```

All four tests, the demo, and both examples passed again.

No fast-math, finite-math-only, or equivalent floating-point assumptions were used.

## Deterministic test scope

- `test_extmat`: callback construction, type check, dimensions, forward/transpose products,
  dense materialization, and column-major vectorization.
- `test_propack`: fixed diagonal singular values for dense and operator paths plus singular
  vector orthogonality.
- `test_trlan`: fixed left-normal singular values, fixed symmetric eigenvalues, dense and
  matrix-free paths, and warm-start validation/plumbing.
- `test_complex`: fixed complex diagonal singular values and Hermitian left-vector
  orthogonality.

## Source and release audit

`scripts/audit.py` passed and checks:

- valid translation counts and matching README/manifest coverage;
- ASCII maintained source and free-form line length;
- `implicit none` in every maintained Fortran unit;
- one documented declaration per dummy argument with explicit `INTENT` or `VALUE`;
- no semicolon-separated executable statements;
- no `double precision`, `real*8`, `d0`/`D0`, or other disallowed maintained real kinds;
- no self-comparison NaN idioms;
- no duplicate maintained Fortran source files;
- no copied RSpectra, rfortran-linalg, ARPACK, BLAS, LAPACK, PROPACK, or nuTRLan sources;
- no objects, modules, executables, caches, nested ZIPs, or build trees in the release.

The maintained package uses one `dp = real64` definition in `src/svd_kinds.f90`.

## Required target-repository FPM validation

With this `svd` directory extracted at the root of `Fortran-from-R-packages`, beside the
existing `RSpectra`, `rfortran-linalg`, and their sibling dependencies, run:

```text
fpm build
fpm test
fpm run svd_demo
fpm run --example dense_svd
fpm run --example matrix_free_svd
fpm clean --all
```

The literal `fpm build`, `fpm test`, and `fpm clean --all` commands could not be executed
in this translation environment because `fpm` was unavailable. A direct attempt to obtain
the current FPM release binary also failed because the compiler container has no external
network access. Consequently, no Windows execution is claimed here; `scripts/validate.bat`
is included for the requested Windows/gfortran target.

## Final archive verification

The release ZIP was created with exactly one top-level `svd/` directory, extracted into
an empty directory, and compared byte-for-byte with the audited source tree. The static
audit passed again on the extraction. The maintained modules, all four tests, the demo,
and both examples were then rebuilt and rerun from the extracted source using the checked
GNU Fortran adapter-validation configuration above; all passed. No forbidden build
artifacts or nested archives were present in the extracted package.
