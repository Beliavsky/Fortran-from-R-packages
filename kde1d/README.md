# kde1d - modern Fortran translation

This directory translates the computational parts of the R package **kde1d 1.2.2** to modern free-form Fortran with FPM. The upstream package is by Thomas Nagler and Thibault Vatter and is distributed under the MIT license.

The implementation covers univariate local-polynomial likelihood kernel density estimation of degree 0, 1, or 2, including automatic plug-in bandwidth selection, case weights, one- and two-sided finite supports, data-adaptive endpoint repair, discrete jittering, zero-inflated hurdle models, density/CDF/quantile evaluation, and simulation.

The upstream implementation accelerates binned Gaussian convolutions with Eigen FFT. This translation evaluates the same finite binned Gaussian derivative convolutions directly. With the fixed 401-knot bandwidth-selection grid this avoids adding an FFT dependency, at the cost of higher asymptotic work during bandwidth selection.

## Dependencies

The package reuses the sibling shared package `rfortran-core` for the common real kind and R-like normal-distribution, ordinary type-7 quantile, and weighted Hyndman-Fan type-7 quantile helpers:

```toml
[dependencies]
rfortran-core = { path = "../rfortran-core" }
```

No BLAS or LAPACK library is required and no translated dependency source is vendored here.

Expected repository layout:

```text
Fortran-from-R-packages/
  kde1d/
  rfortran-core/
```

## Build and test

From the `kde1d` directory:

```text
fpm build
fpm test
fpm run --example basic_kde1d
```

The maintained source is free-form Fortran and uses `dp` from `r_kinds` through `rfortran-core`; `dp` is re-exported by `kde1d_api` for callers.

## Minimal use

```fortran
program example
   use kde1d_api, only : dkde1d, dp, kde1d_fit, kde1d_model
   implicit none

   type(kde1d_model) :: fit
   real(dp) :: x(5)
   integer :: ierr

   x = [-1.0_dp, -0.5_dp, 0.0_dp, 0.5_dp, 1.0_dp]
   call kde1d_fit(x, fit, ierr, bw=0.4_dp, deg=1)
   if (ierr /= 0) error stop "kde1d_fit failed"
   print *, dkde1d([0.0_dp], fit)
end program example
```

`kde1d_fit` accepts optional `xmin`, `xmax`, `mult`, `bw`, `deg`, `weights`, `boundary_repair`, and `grid_size`. Omitting `bw` requests the translated plug-in bandwidth selector. The `type_name` argument accepts continuous, discrete, and zero-inflated aliases analogous to the R API.

For ordered-factor data, pass integer category codes and request `type_name="discrete"`. The generic `equi_jitter` also accepts integer codes directly. R factor labels/classes are intentionally not emulated.

## Translation coverage

**Package status: substantial.** **8 of 8 (100%)** exported computational R functions are mapped. This percentage measures whether computational functions have meaningful Fortran mappings; it does **not** mean full R interface, S3, RNG, or bitwise numerical compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `kde1d` | `kde1d_api` | `kde1d_fit` | substantial |
| `dkde1d` | `kde1d_api` | `dkde1d` | substantial |
| `pkde1d` | `kde1d_api` | `pkde1d` | substantial |
| `qkde1d` | `kde1d_api` | `qkde1d` | substantial |
| `rkde1d` | `kde1d_api` | `rkde1d` | substantial |
| `equi_jitter` | `kde1d_api` | `equi_jitter` | substantial |
| `logLik.kde1d` | `kde1d_api` | `kde1d_loglik` | substantial |
| `summary.kde1d` | `kde1d_api` | `kde1d_summary` | substantial |

Excluded from the denominator are presentation-only S3 methods such as `plot.kde1d`, `lines.kde1d`, `points.kde1d`, and `print.kde1d`.

Material compatibility differences include:

- R ordered-factor objects, factor labels, formula-like metadata, S3 classes, printing, and plotting are not reproduced.
- The model is a Fortran derived type rather than an R list/S3 object.
- Direct finite convolution replaces Eigen FFT in the binned KDE derivative calculations. The estimator equations are preserved, but floating-point reduction order differs.
- `rkde1d` uses the Fortran runtime pseudorandom generator and a one-dimensional Sobol/van-der-Corput-style quasi sequence with optional digital shift. It does not reproduce R `runif` or `randtoolbox::sobol` streams/scrambling exactly.
- R warnings/errors are represented by the explicit `ierr` output of `kde1d_fit` where practical.
- `logLik.kde1d`'s `df` attribute is available as `model%edf`, not as an attached R attribute.

See `COVERAGE.md` and `PROVENANCE.md` for additional detail.
