# mgcv-fortran

A modern Fortran 2018 port of a useful **native numerical core** from the R
package `mgcv` 1.9-4. The project is an FPM package and uses the previously
translated R `splines` package as a local FPM dependency.

This release is not a drop-in, feature-complete replacement for R `mgcv`.
`mgcv` is a large modelling system whose formula parser, S3 classes, sparse
engines, plotting methods, `nlme` integration, and many specialized families
are inseparable from the R runtime. This port translates the central dense
computational workflow and documents omitted subsystems explicitly.

## Implemented numerical API

- Dense penalized IRLS for Gaussian, binomial, Poisson, Gamma,
  inverse-Gaussian, fixed-theta negative-binomial, and fixed-power Tweedie
  models.
- Fixed smoothing parameters, GCV, UBRE, and a documented REML-like dense
  criterion with deterministic coordinate search.
- Cubic regression spline, P-spline, cyclic Fourier spline,
  low-rank thin-plate-style 1D/2D radial bases, random-effect bases, and
  tensor-product design/penalty construction.
- Prediction on the link and response scales, coefficient covariance,
  standard errors, effective degrees of freedom, leverages, and deviance
  residuals.
- Penalized constrained least squares with equality and inequality
  projections, plus monotonicity and convexity constraints.
- Dense equivalents of `XWXd`, `XWyd`, `Xbd`, `diagXVXd`, and `ijXVXd`.
- Matrix roots, symmetric eigensystems, rank estimation, banded and
  tridiagonal Cholesky factors.
- `notExp`, `notLog`, `notExp2`, `notLog2`, normal interval probabilities,
  multivariate normal/t simulation and log densities, Tweedie simulation,
  and weighted chi-square probabilities.
- `gamSim` examples 1 and 2 for normal, Poisson, and binary responses.

See [`docs/API_MAP.md`](docs/API_MAP.md) for exact mapping and
[`docs/PORTING_NOTES.md`](docs/PORTING_NOTES.md) for numerical differences.

## Build with FPM

```text
fpm test
fpm run --example basic_gam
```

The dependency is local:

```toml
[dependencies]
splines = { path = "dependencies/splines" }
```

The bundled dependency manifest also declares `name = "splines"`; FPM
requires this package name to match the dependency-table key exactly.

Standalone GNU Fortran scripts are also supplied:

```text
scripts/run_tests.sh
scripts/run_example.sh
```

On Windows, use `scripts\run_tests.bat` and
`scripts\run_example.bat` from a command prompt containing `gfortran`.

## Minimal example

```fortran
use mgcv

call construct_ps_smooth(x, 12, basis, smooth, status)
design = append_columns(reshape([(1.0_dp, i=1,size(x))], [size(x),1]), basis)

allocate(penalties(13,13,1))
penalties(:,:,1) = embed_penalty(smooth%penalties(:,:,1), 2, 13)

family%id = family_gaussian
call gam_fit(design, y, penalties, model, status, &
             family=family, method=method_gcv)
call model%summary()
```

## License and provenance

The attached `mgcv` source declares `GPL (>= 2)`. The bundled `splines`
translation is also GPL-2.0-or-later. This port is therefore distributed under
GPL-2.0-or-later. Original R/C computational sources are retained under
`original/` for attribution and traceability. No plotting code was translated.
