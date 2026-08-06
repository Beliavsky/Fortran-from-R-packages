# YieldCurve-fortran

A modern Fortran 2018 translation of the computational code in version 5.1 of the R package **YieldCurve** by Sergio Salvino Guirreri.

The port implements Nelson-Siegel and Svensson yield-curve estimation and reconstruction without R, `xts`, or external numerical libraries. It uses native Fortran arrays, bounded golden-section optimization, and a reorthogonalized QR least-squares solver.

## Build

With FPM:

```text
fpm test
fpm run --example basic_fit
```

With GNU Fortran:

```text
./tools/build_gfortran.sh
DEBUG=0 ./tools/build_gfortran.sh
```

## Public API

```fortran
use yieldcurve

call nelson_siegel_fit(rate, maturity, coeff, stat, message)
call ns_rates(coeff, maturity, fitted, stat, message)

call svensson_fit(rate, maturity, coeff, stat, message)
call svensson_rates(coeff, maturity, fitted, "spot", stat, message)
call svensson_rates(coeff, maturity, forward, "forward", stat, message)
```

The fitting and rate-generation procedures are generic over one curve (rank-1 rates/coefficients) and a panel of curves (rank-2 arrays). Missing rate observations represented by IEEE NaNs are omitted during each linear regression.

## Scope

Translated:

- Nelson-Siegel factor loadings, fitting, and rate reconstruction
- Svensson spot/forward loadings, fitting, and rate reconstruction
- R-style candidate grids and bounded one-dimensional maximization
- linear least-squares estimation

Omitted as R-specific infrastructure:

- `xts` date indexes and class restoration
- R data-frame/S3 metadata
- plotting calls that appear only in documentation examples

The original `.rda` datasets are retained under `original/data` for provenance but are not parsed by the Fortran library.

## License

GPL-2.0-or-later, matching the original package's `GPL (>= 2)` declaration. See `LICENSE` and `NOTICE`.
