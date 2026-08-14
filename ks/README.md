# ks-fortran 0.1.0

Modern free-form Fortran computational port of the R package **ks 1.15.3**.

The project is an FPM library and contains no R runtime dependency.  Plotting,
S3 object/display plumbing, data-frame/formula handling and other R UI code are
not translated.  The numerical implementation is organized around explicit
Fortran derived types and array procedures.

## Implemented numerical areas

- Gaussian kernel density estimation in arbitrary dimension, including weights,
  point evaluation, log density, random generation and one-dimensional CDF and
  quantiles.
- Kernel density derivatives of arbitrary order, Gaussian derivative tensors,
  derivative-functionals (`kfe`) and contour/grid helpers.
- Normal-reference, least-squares cross-validation, one-dimensional plug-in and
  smoothed-cross-validation bandwidth selectors; a full-matrix data-driven
  plug-in selector is also provided.
- Exact Gaussian-mixture MISE/AMISE/ISE calculations and MISE/AMISE bandwidth
  optimization, mixture density/simulation/modes/moments, and Student-t mixture
  density/simulation.
- Kernel discriminant analysis with priors/posteriors, classification error and
  confusion matrices.
- Kernel cumulative-distribution estimation (`kcde`) in one or more dimensions.
  Multivariate Gaussian rectangle probabilities use the vendored
  `mvtnorm-fortran` implementation.
- Mean-shift clustering and pointwise density-ridge iteration.
- One- and two-dimensional histogram density estimation.
- Boundary-corrected beta-kernel KDE and empirical copula density estimation.
- Deconvolution KDE weights and regularization cross-validation.
- Curvature and significant-modal-region statistics.
- Balloon variable-bandwidth KDE in two dimensions.
- Convex-hull / high-density support utilities.
- Global two-sample KDE testing.
- Generic linear binning and interpolation plus symmetric convolution.
- `vec`, `vech`, inverse mappings, Kronecker powers, symmetrizers, scaling and
  sphering utilities.

See `API_MAPPING.md` for the R-to-Fortran correspondence and deliberate v0.1.0
boundaries.

## Build

The project requires BLAS and LAPACK.  With FPM:

```text
fpm build
fpm test
fpm run --example basic_kde
```

The manifest deliberately uses strict modern language settings:

```toml
[fortran]
implicit-typing = false
implicit-external = false
source-form = "free"
```

All translated source is `.f90`; explicit interfaces are provided for every
BLAS/LAPACK routine called by the library.

## Licensing

Upstream `ks` is offered under `GPL-2 | GPL-3`.  This distribution elects the
GPL-2 option because it vendors `mvtnorm-fortran`, whose supplied translation is
GPL-2.0-only.  The upstream GPL-2 and GPL-3 texts are retained in `licenses/`,
and the vendored dependency retains its own `LICENSE`, `NOTICE`, and provenance
files.

The original package metadata is retained under `upstream/`.
