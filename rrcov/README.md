# rrcov modern Fortran

**Official CRAN title:** Scalable Robust Estimators with High Breakdown Point

A modern Fortran 2018/FPM translation of the computational parts of the R
package `rrcov` 1.7-8, "Scalable Robust Estimators with High Breakdown Point."
The original package is by Valentin Todorov and contributors.

The port focuses on numerical algorithms. R formula handling, S3/S4 classes,
printing, lattice/base graphics, interactive identification, and R runtime
integration are intentionally omitted.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --target demo_rrcov
```

With GNU Fortran only:

```text
./scripts/build_checked.sh
./scripts/build_optimized.sh
```

On Windows with gfortran, run `scripts\build_checked.bat`.

No external compiled numerical library is required. Symmetric eigendecomposition,
positive-definite repair, inversion, and robust scales are implemented in the
package. Normal, chi-square, F, and incomplete beta/gamma calculations
delegate to the shared `rfortran-core` FPM dependency.

## Basic use

```fortran
use rrcov, only : dp, covariance_result, cov_mcd

real(dp) :: x(100, 4)
type(covariance_result) :: estimate

! Fill x(row, variable).
call cov_mcd(x, estimate, alpha=0.5_dp, nsamp=500, seed=123)

print *, estimate%center
print *, estimate%covariance
print *, estimate%distances
```

A common dispatcher is also available:

```fortran
call robust_covariance(x, "ogk", estimate)
```

Recognized covariance method strings are `classic`, `mcd`, `mve`, `ogk`,
`mest`, `sest`, `mmest`, `sde`, and `mrcd`.

## Main translated areas

- Classical covariance and Mahalanobis distances
- MCD, MVE, OGK, multivariate M, S, MM, Stahel-Donoho, and MRCD estimators
- Robust univariate MAD, Qn, Sn, tau scale, medcouple, and adjusted outlyingness
- Classical and robust covariance-based PCA
- Locantore, projection-pursuit/grid, and Hubert-style robust PCA
- Classical and robust LDA and QDA, Linda-compatible and LdaPP-compatible entry points
- One- and two-sample Hotelling T2 tests
- Wilks one-way MANOVA with Bartlett and Rao approximations
- Confusion matrices, ILR transformation, matrix square roots, and correlations

See `API_MAP.md`, `TRANSLATION_COVERAGE.md`, and `PORTING_NOTES.md` for detailed
mapping and differences from R.

## Array conventions

Observations are rows and variables are columns: `x(n, p)`. Group labels are
integer arrays of length `n`. All public real-valued APIs use `dp`, defined as
`kind(1.0d0)`.

## License

GPL-3.0-or-later, matching the upstream package. The unmodified upstream source
snapshot is retained under `upstream/rrcov-master` for license and provenance.
