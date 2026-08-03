# MTS modern Fortran

A modern Fortran translation of the computational core of the R package
`MTS` 1.2.1 for multivariate time-series analysis.

The project uses Fortran Package Manager (FPM) layout and has no required
external libraries. It retains the upstream Artistic License 2.0 and includes
the complete supplied upstream source under `original/MTS-1.2.1/`.

## Implemented areas

- unrestricted, sparse, and refined VAR models
- VAR order selection, forecasts, impulse responses, FEVD, and Granger tests
- VARMA and VMA simulation, conditional Gaussian estimation, forecasts,
  impulse responses, covariance matrices, psi weights, and differencing
- VARX models, multivariate regression, recursive least squares, and
  regression with vector autoregressive errors
- cross-correlation matrices, multivariate portmanteau tests, multivariate
  ARCH tests, rank-based ARCH diagnostics, and standardized-residual tests
- EWMA covariance, DCC, BEKK(1,1), modified-Cholesky volatility, grouped constrained correlations, and common
  volatility components
- principal components, asymptotic PCA, constrained factor models,
  Stock-Watson factor forecasting, and Bayesian VAR estimation
- known-cointegration and Johansen-style VECM estimation
- full and partial missing-observation estimation
- matrix square roots, Kronecker products, vectorization/half-vectorization,
  matrix polynomial products, pi weights, matrix filters, Corner tables,
  extended cross-correlations, Kronecker indicator specifications, and approximate
  Kronecker-index identification
- self-contained random-number, probability, linear-algebra, and optimization
  routines

See [API_MAP.md](API_MAP.md) for the correspondence with upstream R functions
and for specialized routines that are not direct one-for-one ports.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example var_analysis
fpm run --example varma_analysis
fpm run --example volatility_analysis
fpm run --example factor_vecm_analysis
fpm run --example structural_tools
```

The package version is `1.2.1`, which is valid FPM semantic-version syntax.

## Build without FPM

On a POSIX shell with GNU Fortran:

```text
./scripts/validate.sh
./scripts/validate_optimized.sh
```

On Windows, FPM is the simplest build route. The source uses standard free-form
Fortran 2018 and does not depend on POSIX APIs.

## Basic use

```fortran
program example
   use mts
   implicit none
   integer, parameter :: n = 500, k = 2
   real(dp) :: x(n,k), mu(k), phi(k,k,1), sigma(k,k)
   real(dp), allocatable :: forecast(:,:)
   type(var_model) :: fit

   call set_random_seed(1234)
   mu = [0.05_dp,-0.02_dp]
   phi(:,:,1) = reshape([0.5_dp,0.0_dp,-0.1_dp,0.3_dp],[k,k])
   sigma = 0.1_dp*eye(k)

   call simulate_var(n,mu,phi,sigma,x)
   call fit_var(x,1,fit)
   call predict_var(fit,x,3,forecast)
end program example
```

## Design differences from R

The API is matrix-first and uses derived-type results rather than R lists or S3
objects. Dates, data frames, formula parsing, printing methods, interactive
model selection, plotting, and bundled `.rda` data access are intentionally not
part of the Fortran library. Callers provide numerical arrays and integer lag
or period mappings directly.

The attached translations of `fGarch`, `fBasics`, and `mvtnorm` were used as
layout and numerical-convention references. Their GPL-licensed implementation
source is not copied or linked. The numerical kernels needed here are local to
this Artistic-2.0 project.

## License

Artistic License 2.0. See [LICENSE](LICENSE), [NOTICE](NOTICE), and
[ORIGIN.md](ORIGIN.md).
