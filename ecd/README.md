# ecd-fortran

Modern Fortran translation of the computational algorithms in the R package
`ecd` 0.9.2.4, "Elliptic Lambda Distribution and Option Pricing Model".

The project provides a self-contained double-precision numerical library for
elliptic/cusp distributions, elliptic-lambda and SGED models, stable-count and
Lihn-Laplace process distributions, lambda option pricing, LAMP simulation,
statistical fitting, and supporting numerical utilities.

## Included numerical areas

- General elliptic/cusp distributions and lambda-only distributions
- Density, CDF, survival, quantile, random generation, moments, and statistics
- Cubic, trigonometric, symmetric, asymmetric, and general lambda root solvers
- ECLD and SGED constants, moments, incomplete moments, MGFs, IMGFs, and OGFs
- Closed-form lambda-4/quartic MGF, IMGF, OGF, and smile transformations
- O, V, Q, Qp, U-lag, fixed-point, and quartic option operators
- Laplace, standardized Lihn-Laplace, stable-count, SLD, and QSLD models
- LAMP tau generation, stable random walks, and iterative process simulation
- Levy lambda/skewed densities and cumulant/moment transformations
- Black-Scholes prices, implied volatility, and split polynomial option fitting
- Bounded Nelder-Mead estimation for ECD, ECLD, and SLD/QSLD models
- Adaptive quadrature, root solving, incomplete gamma functions, Bessel K,
  Dawson, erfi, erfcx, scaled quartic error functions, and asymptotic 2F0
- Return differences, lag statistics, sample moments, empirical quantiles,
  quantile binning, and tail partitioning

## Deliberately omitted infrastructure

The following parts of the R package are not numerical algorithms and are not
reimplemented as Fortran library procedures:

- S4/S3 object dispatch and formatted printing
- Rmpfr arbitrary-precision object plumbing
- SQLite database management (`ecdb`)
- YAML symbol configuration and bundled market-data loaders
- xts/zoo/date metadata and data-frame enrichment
- plotting, report tables, parallel R orchestration, and plot-only wrappers

The numeric coordinates and calculations underlying applicable diagnostics are
provided through the Fortran routines. The complete original package tree is
retained under `original/ecd-master` for provenance.

## Precision and dependencies

The Fortran library is self-contained and uses `real(dp)`, where
`dp = kind(1.0d0)`. It does not require R, Rmpfr, GSL, stabledist, optimx,
SQLite, BLAS, or LAPACK.

The original package can dispatch selected calculations to arbitrary precision.
This port uses stable double-precision algorithms and explicit status codes, so
extreme parameter combinations can differ from high-precision Rmpfr results.
See `PORTING.md`.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example option_pricing
fpm run --example lamp_simulation
```

## Build with GNU Fortran

Linux/macOS:

```text
./scripts/run_gfortran_tests.sh strict
./scripts/run_gfortran_tests.sh release
```

Windows command prompt:

```text
scripts\run_gfortran_tests.bat strict
scripts\run_gfortran_tests.bat release
```

## Main module

Most users can import the umbrella module:

```fortran
use ecd_api
```

Individual modules can be imported to reduce namespace size. See `API.md`.

## License

The original Artistic License 2.0 is preserved. See `LICENSE`, `NOTICE.md`, and
the retained original source.
