# tmvtnorm-fortran

Modern free-form Fortran/FPM translation of the computational core of the R package **tmvtnorm 1.7** by Stefan Wilhelm and Manjunath B G.

The project ports the numerical/statistical algorithms and intentionally omits R-specific S3/S4 objects, formula/model-frame machinery, plotting, demos, and dynamic-registration glue.

## Implemented computational API

- Truncated multivariate normal joint density and rectangle probability
- One-dimensional marginal density/CDF/quantile
- Bivariate marginal density
- First and second moments of the doubly truncated multivariate normal
- Johnson-Kotz reduction when only a subset of coordinates is truncated
- Normal rejection sampling
- Dense covariance and precision Gibbs sampling
- General linear-constraint Gibbs sampling for `a <= D x <= b`, including `r > d`
- Sparse precision Gibbs sampling in CSC and triplet form
- Truncated multivariate Student-t density/probability
- Student-t rejection and Gibbs samplers
- Maximum-likelihood estimation of truncated-normal mean/covariance
- Lee moment conditions
- Manjunath-Wilhelm mean/covariance moment conditions
- Two-step GMM fitting using those moment conditions

The attached `mvtnorm-fortran` translation is vendored as an FPM path dependency and supplies multivariate normal/t probabilities, densities, random generation and core linear algebra. The algorithms are not reimplemented here.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The project links no additional libraries beyond what the vendored `mvtnorm-fortran` package itself requires.

## Main module

```fortran
use tmvtnorm
```

Example:

```fortran
real(dp) :: mu(2), sigma(2,2), lower(2), upper(2)
real(dp), allocatable :: x(:,:)
type(tmvnorm_moments_t) :: mom

mu = 0.0_dp
sigma = reshape([1.0_dp,0.6_dp,0.6_dp,1.5_dp],[2,2])
lower = [-1.0_dp,-0.5_dp]
upper = [ 1.5_dp, 2.0_dp]

call mtmvnorm(mu,sigma,lower,upper,mom)
x = rtmvnorm(1000,mu,sigma,lower,upper,algorithm=algorithm_gibbs,burnin=100,seed=1234)
```

## Numerical design choices

The port uses mathematically equivalent block/precision formulations where they reduce duplicated code:

- covariance-form Gaussian Gibbs sampling uses the precision conditional identity after one covariance inversion;
- one- and two-dimensional marginal densities use exact conditional-normal block distributions instead of manually reconstructing the equivalent partial-correlation formulas;
- constrained covariance estimation can use a log-diagonal Cholesky parameterization to keep every optimizer iterate positive definite.

These choices are documented in `PORTING_NOTES.md`.

## Licensing

- Code translated from **tmvtnorm**: GPL-2.0-or-later, matching upstream `GPL (>= 2)`.
- Vendored **mvtnorm-fortran** dependency: GPL-2.0-only under its existing license and attribution.

See `LICENSE`, `LICENSES/`, `NOTICE.md`, and the vendored dependency's own license files.
