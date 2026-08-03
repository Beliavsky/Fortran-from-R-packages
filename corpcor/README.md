# corpcor-fortran

Modern Fortran 2018 translation of the computational API of the R package
`corpcor` 1.6.10.

The package implements weighted moments, covariance/correlation shrinkage,
precision and partial-correlation estimation, compact singular-value
decomposition, pseudoinverses, symmetric matrix powers, rank/condition
checks, positive-definite repair, and symmetric-matrix conversion tools.
It is self-contained and has no external BLAS or LAPACK requirement.

## FPM

```text
fpm build
fpm test
fpm run
```

The public facade is the module `corpcor`. Upstream R names are available with
periods replaced by underscores, for example `cov_shrink`, `pcor_shrink`,
`cor2pcor`, `mpower`, and `wt_moments`. More descriptive names are exported as
well.

Most shrinkage functions return a typed result containing the numerical value,
the shrinkage intensities, flags indicating whether each intensity was
estimated, and a status code.

## Direct compiler scripts

Unix-like systems:

```text
./scripts/build_checked.sh
./scripts/build_optimized.sh
```

Windows with GNU Fortran:

```text
scripts\build_checked.bat
scripts\build_optimized.bat
```

## Matrix orientation

Data matrices use observations by variables, matching the upstream R package.
Fortran matrices are column-major, but all public routines use conventional
mathematical row/column indexing.

## License

GPL-3.0-or-later. The complete retained upstream source is under
`upstream/corpcor-master`.
