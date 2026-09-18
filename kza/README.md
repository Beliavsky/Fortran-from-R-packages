# kza - modern Fortran translation

This directory translates the computational code of the R package **kza 4.2.0** (Kolmogorov-Zurbenko Adaptive Filters) to modern free-form Fortran with FPM packaging. The numerical implementation covers KZ/KZA filtering in one through three dimensions, KZA local variance, KZ Fourier transforms and related spectral calculations, and rolling local variance.

The translation preserves the upstream GPL licensing and provenance. Plotting, printing, R time-series/class construction, interactive output, and other presentation-only R behavior are intentionally omitted.

## Build

When this directory is placed at the root level of `Fortran-from-R-packages`, build and test with:

```text
fpm build
fpm test
fpm run --example basic_kza
```

The translated package has no BLAS, LAPACK, FFTW, or other external numerical-library dependency. In particular, the exported `periodogram` routine uses a dependency-free direct DFT with R's unnormalized forward-transform convention rather than calling FFTW or `stats::fft`.

## Public API

Import the package API with:

```fortran
use kza_api, only : dp, kz, kza, kzsv, kzs, kzft, kztp, periodogram, transfer_function, rlv
```

`kz`, `kza`, and `rlv` are generic interfaces with rank-specific implementations for 1D, 2D, and 3D real arrays. The public KZ/KZA `m` argument is the full R-level window width: internally the per-side radius is `floor(m/2)`, matching kza 4.2.0. Optional `m2`/`m3` arguments expose anisotropic 2D/3D windows in an idiomatic Fortran form.

KZA reproduces the computational options added or corrected in upstream 4.2.0: `min_size`, derivative tolerance, maximum or 99th-percentile adaptive normalization with zero-quantile fallback, 1D tail NA marking, iterative feeding of each KZA pass into the next, and optional four-rotation symmetrization for matrices. Fortran returns the numerical filtered arrays rather than the heterogeneous R S3 list objects.

`kzft` follows upstream's missing-aware local Fourier averaging, including the asymmetric head/tail construction for even or non-integer window values. `periodogram` returns an `n/2 x 2` real matrix containing the same x-axis convention and Fourier magnitudes as the R function. `kztp` returns the complex third-order periodogram matrix.

`rlv` implements both upstream boundary policies. `pad="clamp"` uses only cells actually present at an edge; `pad="zero"` reproduces conceptual zero padding. `krnl=1` preserves the legacy one-sided maximum-over-corners definition.

## Translation coverage

Package status: **substantial**.

Coverage is **9 of 9 (100%)** exported computational R functions. This fraction measures functions with meaningful Fortran mappings; it does **not** mean complete compatibility with R objects, S3 dispatch, attributes, printing, plotting, exact error text, or identical floating-point reduction order.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `kz` | `kza_filters` | `kz` | substantial |
| `kza` | `kza_filters` | `kza` | substantial |
| `kzsv` | `kza_filters` | `kzsv` | substantial |
| `kzs` | `kza_spectral` | `kzs` | substantial |
| `transfer_function` | `kza_spectral` | `transfer_function` | complete |
| `kzft` | `kza_spectral` | `kzft` | substantial |
| `kztp` | `kza_spectral` | `kztp` | substantial |
| `periodogram` | `kza_spectral` | `periodogram` | substantial |
| `rlv` | `kza_rlv` | `rlv` | complete |

The registered `summary.kzp` and `summary.kzsv` methods are presentation/reporting methods and are excluded from the computational coverage denominator under the repository's coverage convention. Plot methods are likewise excluded. Internal, unexported experimental helpers such as `kzp`, `smooth.kzp`, `variation.kzp`, `nonlinearity.kzp`, and `coeff` are not counted as exported computational API mappings.

See `COVERAGE.md` for compatibility notes and `VALIDATION.md` for the exact validation performed in the translation environment.

## Numerical compatibility notes

- KZ moving averages omit non-finite values as the upstream compiled implementation does.
- KZA difference metrics and head/tail allocation are computed once from the KZ baseline, while each adaptive iteration consumes the previous iteration's output, matching the corrected upstream 4.2.0 C code.
- Dimension-one 2D/3D edge cases use the intended zero directional derivative rather than reproducing undefined out-of-bounds behavior.
- The direct DFT in `periodogram` is mathematically equivalent to the R FFT convention used by the package, but is `O(n^2)` and can differ in the last floating-point bits because the summation order differs.
- `kzsv` accepts the numeric KZA series, KZ baseline, and window controls directly because there is no R S3 object in the Fortran API.
- R `ts` metadata, calls, list components, plotting, console progress messages, and summary printing are outside the translated computational surface.

## Licensing and provenance

The upstream package declares GPL-3 licensing. Original computational R and C sources used for the translation are retained under `upstream/` with their original notices and comments. See `NOTICE.md` and `PROVENANCE.md`.
