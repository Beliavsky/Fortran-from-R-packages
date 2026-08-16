# DiscreteLaplace-fortran

Modern Fortran 2018/FPM translation of the computational code in the R package
`DiscreteLaplace` 1.1.1 by Alessandro Barbiero and Riccardo Inchingolo.

The port is standalone: base-R optimization has been replaced by a small
Nelder-Mead implementation.  There are no external runtime dependencies.

## Main API

The original exported names are retained where Fortran permits them:
`ddlaplace`, `pdlaplace`, `qdlaplace`, `rdlaplace`, `ddlaplace2`,
`palaplace2`, `pdlaplace2`, `qdlaplace2`, `rdlaplace2`, `edlaplace`,
`edlaplace2`, `ifi`, `ifi2`, `iofi2`, `estdlaplace`, `estdlaplace2`,
`dlaplacelike2`, and `loss`.

`edlaplace`, `edlaplace2`, and `estdlaplace` return derived types because the
R functions return lists.  Quantiles/random draws return `integer(int64)`.

## Build

```text
fpm build
fpm test
```

The release was also checked directly with gfortran 14.2 using Fortran 2018,
`-Werror=implicit-interface`, and `-fcheck=all`.

See `API_MAP.md` and `PORTING_NOTES.md` for interface and compatibility notes.
