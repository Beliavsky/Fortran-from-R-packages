# Validation

Validation was performed on 2026-09-20 with GNU Fortran 14.2.0.

## Dependency review

Before implementing shared functionality, the current `Fortran-from-R-packages` repository was checked. The translation reuses these existing sibling packages:

- `rfortran-core` (`name = "rfortran-core"`) for `r_kinds::dp`, normal and chi-square distribution functions, average ranks, and sample variance.
- `SuppDists` (`name = "suppdists-fortran"`) for `suppdists::norm_order`.

Their current public procedure interfaces were checked against the calls in this package. No dependency source is copied into `kSamples`.

## Compiler validation

The package sources, tests, and example were compiled twice in dependency order against lightweight API-compatible dependency stubs kept outside this package. The stubs were used only because the sandbox does not contain the sibling repository directories; they are not included in the release ZIP.

The runtime-checking build used:

```text
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace -Wall -Wextra -Werror=line-truncation
```

The optimized build used:

```text
gfortran -std=f2018 -O2 -Werror=line-truncation
```

Both builds produced the same successful test/example results:

```text
basic kSamples tests passed
exact kSamples tests passed
QN statistic = 3.5769231, exact p = 0.1729870
JT statistic = 9.0000000, exact p = 0.1666667
```

The deterministic tests cover convolution, exact JT distribution functions, QN, both Anderson-Darling statistic versions, JT, Steel exact allocation, contingency tables, blocked combinations, and the asymptotic Steel confidence-interval path. Simulation tests are not used for release pass/fail so validation does not depend on RNG state.

## FPM limitation in this sandbox

The required commands were attempted literally:

```text
fpm build
fpm test
fpm run --example rank-tests
fpm clean --all
```

All four returned exit status 127 because `fpm` is not installed in the sandbox (`fpm: not found`). An attempt to install the official PyPI `fpm` package also failed because command-line network/DNS access is disabled. Therefore this document does not claim that FPM itself was executed successfully here.

The `fpm.toml` manifest was parsed independently as TOML, its path dependencies were checked against the current sibling package names, and the same package sources/tests/examples were compiled and run directly with gfortran as described above.

## Release audits

Before packaging, automated static checks verified:

- 15 regular exports in `upstream/NAMESPACE`, of which `pp.kSamples` is plotting-only and excluded; all 14 exported computational functions have exactly one coverage mapping.
- `README.md` and `fpm.toml` both report 14 of 14 (100%) mapped computational functions and package status `substantial`.
- Every Fortran dummy argument has `INTENT` or `VALUE`, is declared on its own line, and has a meaningful trailing FORD `!!` comment.
- All maintained Fortran is free-form ASCII and uses the single `dp` kind imported from `r_kinds`.
- No `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals occur in maintained Fortran.
- No disallowed semicolon-separated Fortran statements occur.
- No self-comparison NaN tests occur in this package.
- No system BLAS/LAPACK link flags occur.
- No copied `r.f90`, `r_mod.f90`, BLAS, LAPACK, ARPACK, `rfortran-core`, or `SuppDists` implementation source occurs.
- No duplicate maintained Fortran source files occur.
- No object files, module files, executables, caches, build directories, or ZIP files are included in the package tree.
