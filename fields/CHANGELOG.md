# Changelog

## 0.1.0

Initial computational translation of fields 17.3.

- Added modern distance, covariance, polynomial, spline, Kriging, thin-plate spline, sparse Kriging, fast Tps, spatial-process fitting, statistics/variogram, grid, FFT/circulant, geometry, simulation, and quantile-smoothing modules.
- Converted the upstream fixed-form `fieldsF77Code.f` numerical kernels to free-form module source and wrapped them with typed APIs.
- Vendored the user-supplied `spam-fortran-v0.1.0` as an FPM path dependency.
- Added a local spam identity-pivot initialization fix required for `pivot='none'` sparse Cholesky.
- Added explicit initialization of native spline error flags before every `css`/`rcss` call.
- Added 10 regression programs and a basic example.
- Retained upstream fields source/metadata and dependency notices for license/provenance review.
