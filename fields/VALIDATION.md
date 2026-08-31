# Validation

Validation was performed in the build environment with **GNU Fortran 14.2.0**, system BLAS, and LAPACK.

FPM itself was not installed, so the FPM source/dependency graph was compiled directly. The root `fpm.toml` parses successfully as TOML and retains the bundled `spam-fortran` path dependency.

## Compiler settings

Modern `fields` modules and test/example programs were compiled with the equivalent of:

```text
-O0 -g -std=f2018 -ffree-line-length-none -fcheck=all -Werror=implicit-interface
```

The mechanically converted legacy native source and inherited spam native sources were compiled in free form with runtime checking where applicable but without forcing all inherited intra-file calls to acquire new explicit interfaces. Linking used `-llapack -lblas`.

## Regression programs

All of the following passed in a clean rebuild from the packaged source tree:

- `test_additional.f90` — Matérn range conversion, Kriging derivative, weighted/group one-way utilities.
- `test_covariance_distance.f90` — covariance and distance kernels.
- `test_fft_interp_wendland.f90` — Fourier interpolation at source knots and native Wendland-grid multiplication.
- `test_fft_spatial_process.f90` — circulant embedding and covariance-parameter fitting.
- `test_grid_stats_variogram.f90` — grid interpolation, summaries/binning, and variograms.
- `test_kriging.f90` — dense universal Kriging and prediction.
- `test_quantile.f90` — `QSreg`/`QTps`-style quantile smoothers and loss functions.
- `test_sparse.f90` — sparse spam-backed Wendland Kriging versus dense reference behavior, including the no-pivot compatibility path.
- `test_spline.f90` — cubic smoothing spline interpolation/prediction/derivatives and df-to-lambda inversion.
- `test_tps.f90` — thin-plate spline fitting/prediction.

`example/basic.f90` also built and ran successfully.

## Important validation fixes

Two defects were found by the regression pass and corrected:

1. The vendored spam `pivot='none'` wrapper failed to initialize its identity permutation, causing a segmentation fault inside the inherited sparse Cholesky routine.
2. The native `css` routine does not assign `ierr=0` on a successful return. Modern wrappers now initialize the flag before every `css`/`rcss` call; otherwise a nonlinear prediction could inspect indeterminate stack data as an error code.
