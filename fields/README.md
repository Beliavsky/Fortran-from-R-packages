# fields-fortran

`fields-fortran` is a modern free-format Fortran/FPM translation of the computational core of the R package **fields 17.3** (2026-05-03), with the user-supplied `spam-fortran` port bundled as an FPM path dependency.

The port is aimed at numerical use from Fortran. It translates or wraps the package's spatial-statistics algorithms rather than reproducing R's S3/formula/graphics layer.

## Implemented numerical areas

- Euclidean and great-circle distances, compact-neighbor conversions, diagonal updates.
- Exponential, Gaussian, powered-exponential, Matérn, cubic, radial/thin-plate, Wendland/taper, and Paciorek covariance calculations.
- Matérn correlation-to-range conversion and anisotropic stationary covariance matrices.
- Polynomial/null-space construction and polynomial gradients.
- Cubic smoothing splines, GCV selection, df-to-lambda inversion, derivatives, and the original robust `rcss` spline kernel.
- `QSreg`- and `QTps`-style iterative quantile smoothing with asymmetric pseudo-responses and pseudo-data leave-one-out criteria.
- Dense universal Kriging with GCV, ML and REML profiles, prediction covariance/standard errors, derivative prediction, conditional simulation, and multiple responses.
- Thin-plate splines with exact radial/null-space saddle systems and GCV.
- Sparse `mKrig`/Wendland Kriging and `fastTps` using the bundled `spam` sparse Cholesky/solve implementation.
- `spatialProcess`-style covariance parameter profiling and bounded joint optimization of range, nugget ratio, and optional Matérn smoothness.
- Ordinary and Cressie robust variograms, descriptive/bin statistics, and weighted one-way group summaries.
- Bilinear interpolation, grid boxes, off-grid covariance interpolation weights, direct and FFT image smoothing, Fourier surface interpolation, and native Wendland-grid multiplication.
- 2-D circulant-embedding Gaussian-field simulation, general dense random-field simulation, Paciorek-field simulation, spatial-data simulation, and conditional Kriging simulation.
- Polygon membership and minimax/greedy space-filling utilities.

The original `fieldsF77Code.f` numerical routines are retained in `upstream/` and mechanically converted to free-form Fortran in `src/native/fields_native.f90`; typed modern modules are layered over them.

## Build

The project uses BLAS and LAPACK and has the bundled spam port as a local FPM dependency:

```text
fpm build
fpm test
```

FPM was not available in the validation container. The same source graph was therefore compiled directly with GNU Fortran 14.2.0, BLAS and LAPACK. See `VALIDATION.md`.

The public convenience module is:

```fortran
use fields
```

A small example is in `example/basic.f90`.

## Deliberately not translated

Plotting, palettes, maps, R S3 methods, formula/model-frame construction, printing/summary formatting, datasets, and other R-specific presentation or dispatch code are excluded. Some R performance orchestration (for example, object-caching wrappers around fast grid prediction) is represented by lower-level reusable numerical kernels instead of R-shaped objects.

See `API_MAPPING.md` and `PORTING_NOTES.md` for detailed coverage and design choices.

## License and provenance

The translated `fields` code is GPL-2.0-or-later, following the upstream package declaration. The bundled `spam` dependency has multiple upstream provenance/licensing layers and must **not** be treated as uniformly GPL. See `NOTICE.md`, `upstream/LICENSE.note`, and `vendor/spam/NOTICE.md` before redistribution, especially for commercial use.
