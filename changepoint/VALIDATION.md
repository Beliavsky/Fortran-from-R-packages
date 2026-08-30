# Validation record

Validation was performed in the translation sandbox on 2026-08-29.

## Toolchain

- GNU Fortran: `GNU Fortran (Debian 14.2.0-19) 14.2.0`
- LAPACK/BLAS: system libraries used for the shared `rfortran-linalg` dependency
- FPM: not installed in this sandbox

## Shared dependency validation

The package is configured only with sibling path dependencies:

```toml
rfortran-core = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

The translation was compiled and linked against the actual shared
`rfortran-core/src/r_kinds.f90` and `rfortran-linalg/src/r_linalg.f90` sources
available in the validation workspace. The shared `r_linalg` SVD least-squares
routine linked to system LAPACK/BLAS. No dependency source is included in this
release directory.

## Strict compiler validation

Every maintained `changepoint` library source, the deterministic test program,
and both examples were compiled with the package sources under:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Dependency modules were compiled separately as dependencies. The package build
was warning-free under the strict flags above.

Test output:

```text
all changepoint tests passed
```

Example output:

```text
number of changepoints: 2
locations: 4 8
number of changepoints: 1
locations: 6
```

The deterministic suite covers likelihood AMOC, genuinely pruned PELT, Binary
Segmentation, exact Segment Neighborhood, PELT-versus-exact-DP agreement,
Normal/Exponential/Gamma/Poisson costs, zero-valued Exponential input,
noninteger-Poisson rejection, CUSUM/CSS methods, CROPS, regression AMOC/PELT,
segment parameter fits, penalty calculation, the developer decision helper,
and asymptotic AMOC transforms.

## FPM gate

The required FPM commands were attempted from the package root immediately
before release preparation. The sandbox does not provide an `fpm` executable:

```text
+ fpm build
bash: line 1: fpm: command not found
exit=127

+ fpm test
bash: line 1: fpm: command not found
exit=127

+ fpm clean --all
bash: line 1: fpm: command not found
exit=127
```

Therefore this record does **not** claim that FPM itself was executed
successfully in this environment. `fpm.toml` was parsed as TOML and its sibling
dependency paths and free-form Fortran settings were checked. On a checkout
that has FPM installed, `python tools/release_check.py` executes `fpm build`,
`fpm test`, all release-tree scans, and `fpm clean --all`.

## Release-tree hygiene

The final pre-archive scans reported:

- 15 maintained free-form `.f90` files
- 15 unique SHA-256 hashes for those 15 files
- maximum maintained Fortran line length: 121 columns
- every maintained Fortran file imports the same `dp` from `r_kinds`
- the public `changepoint` module re-exports `dp`
- no `double precision`, `real*8`, `real64`, package-local `iso_fortran_env`,
  `kind(0.0d0)`, or D-exponent real literals
- no semicolon-separated Fortran statements
- no vendored `rfortran-core`, `rfortran-linalg`, BLAS, LAPACK, ARPACK,
  `r.f90`, `r_mod.f90`, `r_kinds.f90`, or `r_linalg.f90`
- no duplicate maintained Fortran source files
- no object/module/executable/cache/build/ZIP artifacts in the release tree

`tools/release_check.py` encodes these checks for repeatable release use.
