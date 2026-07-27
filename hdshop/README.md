# hdshop-fortran

A self-contained modern Fortran translation of the numerical algorithms in
HDShOP 0.1.7, packaged for the Fortran Package Manager (FPM).

HDShOP constructs shrinkage estimators for high-dimensional mean-variance and
global-minimum-variance portfolios. The original package is retained under
`original/HDShOP-0.1.7` and remains the authoritative source for the published
methodology and citations.

## Features

- Sample covariance with finite-observation handling.
- Bodnar-Golz-Parolya linear covariance shrinkage.
- Ledoit-Wolf analytical nonlinear covariance shrinkage.
- Bodnar-Gupta-Parolya linear precision-matrix shrinkage.
- Bayes-Stein, James-Stein, and BOP19 mean shrinkage.
- Traditional mean-variance and GMV portfolios.
- BDOPS21 and BDPS19 shrinkage portfolios for both `p < n` and `p > n`.
- Moore-Penrose inverse for singular high-dimensional sample covariances.
- Per-weight asymptotic intervals and tests when `p < n`.
- Mean-variance portfolio-efficiency testing.
- Bayesian efficient-frontier arrays.
- Seeded random covariance matrices with prescribed eigenvalues.
- Upstream-style compatibility entry points and descriptive typed APIs.

The library has no external numerical dependencies.

## Build

```text
fpm build
fpm test
fpm run hdshop_demo
fpm run --example shrinkage_estimators
fpm run --example high_dimensional_portfolio
```

Direct validation is also provided:

```text
scripts/validate.sh
scripts/validate.sh optimized
```

On Windows with GNU Fortran:

```bat
scripts\validate.bat
```

## Basic use

```fortran
use hdshop

real(dp) :: returns(assets, observations)
real(dp) :: target(assets)
type(portfolio_result) :: fit

target = 1.0_dp / real(assets, dp)
fit = shrinkage_gmv_portfolio(returns, target)

if (fit%ok) print *, fit%weights
```

Input returns follow the R package convention: rows are assets and columns are
observations.

## API conventions

R lists and S3 objects are represented by derived types. Plotting is replaced by
`bayesian_frontier`/`frontier_data`, which return the numerical curve and
portfolio points. Fortran uses zero status flags and messages instead of R
exceptions.

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for details.

## License

GPL-3.0-only, matching HDShOP 0.1.7. See `LICENSE` and `NOTICE`.
