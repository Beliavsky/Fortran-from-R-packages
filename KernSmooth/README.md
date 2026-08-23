# KernSmooth-fortran

Modern Fortran/FPM translation of the numerical core of R's `KernSmooth` package.

## Implemented public numerical routines

- `linbin`, `linbin2d`, `rlbin`
- `bkde`, `bkde2d`
- `bkfe`
- `locpoly`
- `sdiag`, `sstdiag`
- `dpih`, `dpik`, `dpill`

The module uses `real(dp)` with `dp = kind(1.0d0)` and free-form Fortran 2018 source.

## Notes

The original package uses FFT convolution for binned density estimates. This port uses direct grid convolution, which is mathematically equivalent but slower for very large grids and avoids an FFT dependency.

The original low-level fixed-form Fortran sources are retained under `orig/` for provenance and algorithm comparison. They are not part of the FPM build.

The `dpill`, `sdiag`, and `sstdiag` implementations are modern direct local-polynomial equivalents of the original blocked/LINPACK machinery rather than line-for-line translations of those internal work arrays.

Build with:

```sh
fpm test
```
