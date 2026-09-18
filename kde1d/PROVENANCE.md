# Provenance and licensing

## Upstream R package

- Package: `kde1d`
- Version in the supplied source archive: 1.2.2
- Authors: Thomas Nagler and Thibault Vatter
- Upstream project URL recorded in `DESCRIPTION`: `https://tnagler.github.io/kde1d/`
- Upstream source repository recorded in `DESCRIPTION`: `https://github.com/tnagler/kde1d`
- License: MIT + upstream R `LICENSE` metadata file

The original R license metadata is preserved verbatim as `upstream/LICENSE`; it is also retained at package root as `LICENSE`. `LICENSE-MIT` supplies the full MIT grant for tooling that does not interpret R's `YEAR` / `COPYRIGHT HOLDER` metadata convention.

The upstream package includes the MIT-licensed `kde1d-cpp` computational implementation. Its license is preserved verbatim at `upstream/include/kde1d-cpp/LICENSE`. The maintained Fortran source identifies the upstream authors and the translated numerical provenance in source headers.

The exact R files used for exported-function coverage are retained under `upstream/R/`, and the supplied `DESCRIPTION` and `NAMESPACE` are retained under `upstream/`.

## Numerical provenance

`src/kde1d_bandwidth.f90`, `src/kde1d_interpolation.f90`, and `src/kde1d_api.f90` translate numerical behavior from the R/C++ implementation in the supplied archive, including its plug-in bandwidth formulas, interpolation scheme, support transformations, local-polynomial corrections, jittering, boundary repair, discrete normalization, and zero-inflated behavior.

The upstream C++ `KdeFFT` uses Eigen's FFT to evaluate binned Gaussian convolutions. The Fortran translation deliberately uses direct finite convolution on the same regular grids and with the same Gaussian derivative kernels/truncation rules. This avoids introducing an FFT dependency and changes floating-point reduction order and computational complexity, but not the intended estimator equations.

## Shared Fortran dependency

The package depends on the existing sibling `rfortran-core` package from `Fortran-from-R-packages` for `dp`, standard-normal density/CDF/quantile helpers, medians, R type-7 quantiles, and the weighted Hyndman-Fan type-7 quantile convention used by kde1d. Dependency source is not copied into this package.

No BLAS, LAPACK, ARPACK, Eigen, Boost, Rcpp, or randtoolbox source is copied or linked by the translated package.
