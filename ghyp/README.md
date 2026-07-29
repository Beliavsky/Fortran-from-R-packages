# ghyp-fortran

A dependency-free modern Fortran/FPM numerical-core translation of the R package
`ghyp` 1.6.5.

The library supports the univariate and multivariate generalized hyperbolic
family and its principal special cases:

- hyperbolic;
- normal inverse Gaussian;
- variance gamma;
- symmetric and skewed Student t;
- Gaussian.

It also includes the generalized inverse Gaussian mixing distribution,
likelihood fitting, risk measures, transformations, ES attribution, and
portfolio optimization.

## Build

```text
fpm build
fpm test
fpm run ghyp_demo
fpm run --example distributions_and_risk
fpm run --example fitting_and_portfolio
```

The public module is:

```fortran
use ghyp
```

All calculations use double precision through:

```fortran
integer, parameter :: dp = kind(1.0d0)
```

## Basic distribution example

```fortran
use ghyp

implicit none

type(ghyp_model_type) :: model
type(moments_result) :: moments

model = ghyp_ad( &
   lambda = 0.7_dp, &
   alpha = 1.8_dp, &
   delta = 1.2_dp, &
   beta = [0.3_dp], &
   mu = [0.2_dp], &
   delta_matrix = reshape([1.0_dp],[1,1]) &
)

if (.not. model%ok) error stop trim(model%message)

moments = ghyp_moments(model)

print *, dghyp(0.5_dp,model)
print *, pghyp(0.5_dp,model)
print *, qghyp(0.95_dp,model)
print *, esghyp(0.95_dp,model,loss=.true.)
print *, moments%mean
print *, moments%covariance
```

## Constructors

The main constructors are:

```fortran
model = ghyp_uv(lambda,chi,psi,mu,sigma,gamma)
model = ghyp_mv(lambda,chi,psi,mu,scatter,gamma)
model = ghyp_ad(lambda,alpha,delta,beta,mu,delta_matrix)

model = hyp_uv(...)
model = nig_uv(...)
model = student_t_uv(...)
model = vg_uv(...)
model = gaussian_uv(...)
model = gaussian_mv(...)
```

For the univariate `chi/psi` constructor, `sigma` is a standard-deviation
scale. The multivariate constructor accepts a positive-definite scatter
matrix.

## GIG distribution

The translated generalized inverse Gaussian API includes:

```fortran
dgig
log_dgig
pgig
qgig
rgig
esgig

gig_raw_moment
gig_mean
gig_variance
gig_mean_log
gig_mean_inverse
```

Gamma and inverse-gamma limiting cases are handled directly.

## GH distribution and risk measures

The main routines include:

```fortran
dghyp
log_dghyp
pghyp
pghyp_rectangle
qghyp
rghyp

ghyp_moment
ghyp_skewness
ghyp_kurtosis
esghyp
ghyp_omega
esghyp_attribution
```

`pghyp_rectangle` follows the upstream package's simulation-based approach for
multivariate probabilities and returns a Monte Carlo standard error.

## Transformations and utilities

```fortran
transformed = transform_ghyp(model,multiplier,summand)
marginal = subset_ghyp(model,[1,3])
standardized = standardize_ghyp(model)
parameters = ghyp_alpha_delta(model)
qq = qqghyp_data(data,model)
```

Linear transformations preserve the normal mean-variance-mixture
representation exactly.

## Fitting

```fortran
fit = fit_ghyp_uv(data,"nig")
fit = fit_ghyp_mv(data,"ghyp")
```

Available family names are:

```text
ghyp, hyp, nig, student, vg, gaussian
```

The fit result contains the fitted model, log likelihood, AIC, BIC,
transformed optimizer parameters, numerical Hessian covariance, standard
errors, iteration count, and convergence status.

Additional routines include:

```fortran
fit_gaussian_uv
fit_gaussian_mv
likelihood_ratio_test
step_aic_ghyp
```

The upstream multivariate implementation uses specialized EM updates. This
port uses a common direct maximum-likelihood objective with native Nelder-Mead
optimization. It targets the same likelihood but is not expected to reproduce
the upstream iteration path or fitted values bit for bit.

## Portfolio optimization

```fortran
portfolio = portfolio_optimize( &
   market_model, &
   risk_measure = "sd", &
   problem = "minimum.risk" &
)
```

Supported risk measures:

```text
sd
value.at.risk
expected.shortfall
```

Supported problems:

```text
minimum.risk
tangency
target.return
```

Mean-variance and symmetric-distribution cases use analytical linear systems.
Non-symmetric VaR and ES cases use native Nelder-Mead optimization.

## Numerical implementation

The library is self-contained. It includes:

- cached Gauss-Legendre quadrature;
- a native real-order modified Bessel K calculation;
- incomplete-beta and Student-t calculations;
- Cholesky and general linear-system solvers;
- deterministic random-number generation;
- Nelder-Mead optimization;
- numerical Hessians.

No R, BLAS, LAPACK, GSL, Boost, SciPy, or external optimization library is
required.

General GIG random generation uses inverse-CDF sampling. It is deterministic
and robust but slower than the specialized rejection sampler in the upstream C
source.

The fitting and non-symmetric portfolio optimizers keep objective context in
module state. Calls to those optimizers should therefore be externally
serialized when used from multiple threads.

## Scope

The numerical portions of the exported package API are represented. R-specific
infrastructure is not compiled, including:

- S4 classes and methods;
- plotting, histograms, pair plots, and graphics devices;
- formatted print and summary methods;
- R data frames and formula handling;
- bundled `.rda` data as compiled inputs;
- exact reproduction of R optimizer control objects.

The original package and data are retained for provenance.

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for details.

## License

`GPL-2.0-or-later`, matching the upstream package.
