# RiskPortfolios Fortran

A modern Fortran 2018 and FPM implementation of the computational methods in
RiskPortfolios 2.1.7.

The library computes mean and semideviation estimates, ten covariance matrix
estimators, and seven risk-based portfolio constructions. It has no R runtime
dependency. BLAS and LAPACK are used for linear solves and symmetric
eigendecompositions.

## Implemented features

Mean estimators:

- Arithmetic sample mean
- Exponentially weighted mean
- Bayes-Stein mean
- Martellini implied mean, represented by sample volatility as in the R source

Covariance estimators:

- Sample covariance
- Exponentially weighted covariance
- Ledoit-Wolf market shrinkage
- Factor covariance
- Constant-correlation covariance
- Ledoit-Wolf constant-correlation shrinkage
- One-parameter shrinkage
- Diagonal shrinkage
- Large-dimensional market shrinkage
- Bayes-Stein covariance

Portfolio methods:

- Mean-variance
- Minimum variance
- Inverse volatility
- Equal risk contribution
- Maximum diversification
- Risk efficient
- Maximum decorrelation

Constraint modes:

- Budget constraint only
- Long only
- User lower and upper bounds
- Gross-exposure limit

## Requirements

- A Fortran 2018 compiler
- Fortran Package Manager (FPM), for the standard build
- BLAS and LAPACK

GNU Fortran example:

```text
fpm build
fpm test
fpm run demo_riskportfolios
```

When BLAS and LAPACK are not found automatically, add the appropriate library
search path for the operating system and compiler.

## Minimal example

```fortran
program example
   use riskportfolios
   implicit none

   real(dp) :: returns(100, 4)
   real(dp), allocatable :: mu(:), sigma(:, :), semidev(:), weights(:)
   type(portfolio_control) :: control
   integer :: info

   ! Fill returns with observations in rows and assets in columns.

   call mean_estimation(returns, mu, MEAN_BAYES_STEIN, info=info)
   call covariance_estimation(returns, sigma, COV_LEDOIT_WOLF, info=info)
   call semideviation_estimation(returns, semidev, SEMIDEV_EWMA, &
      lambda=0.94_dp, info=info)

   control%constraint = CONSTRAINT_LONG_ONLY
   call optimal_portfolio(sigma, weights, PORT_EQUAL_RISK_CONTRIBUTION, &
      mu=mu, semidev=semidev, control=control, info=info)
end program example
```

## Project layout

- `src/`: library source
- `app/`: demonstration executable
- `example/`: additional examples
- `test/`: numerical and behavioral tests
- `licenses/`: complete GPL texts
- `original/`: original package metadata and R source for provenance
- `API.md`: public API reference
- `PORTING.md`: source mapping and documented differences
- `TESTING.md`: validation details

## License and citation

The original package is GPL (>= 2). This port is therefore distributed under
GPL-2.0-or-later. See `LICENSE`, `NOTICE.md`, and `licenses/`.

Publications using this port should retain the original RiskPortfolios citation
described in `NOTICE.md` and `original/CITATION`.
