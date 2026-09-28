# Validation record

Validation was performed on the exact maintained source tree for this checkpoint.

## Compiler validation

Compiler: GNU Fortran 14.2.0.

Two independent builds were performed against temporary, interface-compatible validation shims for the public `rfortran-core` and `rfortran-linalg` APIs. The shims live outside the package tree and are not included in the ZIP. They use small pure-Fortran implementations sufficient to exercise the translated package interfaces in this sandbox; the package itself continues to depend on the repository's sibling `rfortran-core` and `rfortran-linalg` packages and contains no copied BLAS/LAPACK or dependency source.

Strict checked flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -fbacktrace
```

Optimized flags:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface
```

Under both modes all seven deterministic test programs passed:

- `test_statistics`
- `test_covariance`
- `test_factor_simulation`
- `test_models`
- `test_extended`
- `test_remaining`
- `test_remaining_models`

The `example/basic.f90` program also passed in both modes and printed the one-step functional forecast `1.4000 1.8000 3.6000`.

The checked build exercises IEEE-NaN missing-value paths, bounds checking, allocation checks, and the newly translated FAR/VAR, stationarity, density-FPCA, compositional-density, LQD, skew-t, MFDM, and smooth-GEV workflows. Deterministic tests use fixed Fortran RNG seeds where simulation is involved; R and Fortran RNG streams are not claimed to be identical.

## Windows/FPM portability corrections

A downstream Windows `fpm test` run of the preceding checkpoint exposed two portability problems in `test_remaining` and `test_remaining_models`:

- the FTSM PCA used `rfortran-linalg::thin_svd`, reaching the wide-matrix `DGESDD` path and returning a nonzero status in the `farforecast` regression test;
- the transformed skew-t parameter VAR used SVD least squares on a nearly collinear small design, which could return a nonzero status in `skew_t_fun`.

This checkpoint removes both sensitive paths. FTSM PCA now obtains the same weighted principal subspace from the symmetric weighted covariance eigendecomposition. The small score/parameter VAR regression now solves a lightly ridge-stabilized normal equation through `rfortran-linalg::solve_spd`. The skew-t likelihood is evaluated in log-density form so trial simplex points do not create avoidable subnormal densities. With `-ffpe-summary=all`, `test_remaining_models` passes in the sandbox without an IEEE denormal summary.

A subsequent downstream Windows `fpm test` run exposed the same SVD/LAPACK sensitivity in `GAEVforecast`: the deterministic GAEV regression test returned a nonzero status and reported `IEEE_DENORMAL`. The GAEV implementation no longer calls SVD least squares. Its initial spline regression and small trend/VAR regressions now use ridge-stabilized normal equations through `rfortran-linalg::solve_spd`, and extreme GEV exponential/power evaluations use bounded log-domain exponentiation to avoid subnormal underflow. The full checked test suite passes with `-ffpe-summary=all` in the sandbox after this correction.

## Translation coverage

The audited coverage basis contains 39 distinct exported computational R functions. All 39 have meaningful Fortran mappings, so mapped-function coverage is 39 of 39 (100.0%). This fraction measures mapped computational functions, not complete R interface or package-engine compatibility. Several high-level functions are deliberately marked `partial`, notably those whose upstream implementations rely on JAGS/R2jags, `evgam`, `vars`, `forecast`, `fda`, `fda.usc`, or R object/formula machinery.

## FPM availability

The required commands were attempted on the final tree:

```text
fpm build
fpm test
fpm run --example basic
fpm clean --all
```

The execution environment does not have the `fpm` executable installed, so each command returns exit status 127 (`fpm: command not found`). This checkpoint therefore does **not** claim that literal FPM execution occurred in the sandbox. The manifest is structured for the target repository layout with sibling `../rfortran-core` and `../rfortran-linalg` dependencies.

## Release audit

Before packaging, the maintained Fortran source/test/example tree is checked for:

- one dummy argument per declaration with explicit `INTENT` or `VALUE`;
- a trailing `!!` FORD documentation comment on every dummy declaration;
- free-form source and maximum line length of 132 columns;
- ASCII-only maintained Fortran files;
- no semicolon-separated executable statements;
- no `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals;
- no self-comparison NaN tests;
- no unsafe fast-math flags or system BLAS/LAPACK link flags;
- no duplicate Fortran source basenames or duplicate maintained Fortran files;
- no copied dependency source;
- no object, module, executable, archive, cache, or nested ZIP build products.

The `fpm.toml` coverage arithmetic, all mapping entries, and all recorded R/Fortran source paths are also checked for consistency.
