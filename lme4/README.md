# lme4-fortran

A modern Fortran 2018 computational-core port of numerical functionality from the R package **lme4**.

The API accepts numeric response, fixed-effect, random-effect, covariate, and grouping arrays directly. It intentionally does not reproduce R formulas, factors, contrasts, model frames, S3/S4 classes, plotting, or `Matrix`/CHOLMOD objects.

## Implemented

### Linear mixed models

- Gaussian LMMs fitted by ML or REML
- Multiple independent random-effect terms
- Correlated random intercepts and slopes
- Unstructured, diagonal, compound-symmetry, and AR(1) covariance structures
- Dense marginal-covariance solver through `fit_lmm`
- Penalized-least-squares/Woodbury solver through `fit_lmm_pls`
- BLUPs, fixed-effect covariance, fitted values, residuals, log likelihood, deviance, AIC, and BIC

### Generalized mixed models

- Dense PIRLS plus first-order Laplace fitting
- Built-in binomial-logit, Poisson-log, Gamma-log, inverse-Gaussian-log, and negative-binomial-log families
- Fixed dispersion for Gamma and inverse Gaussian
- Fixed or profiled negative-binomial size through `fit_glmer_nb`
- User-defined family/link callbacks through `family_spec_t` and `fit_glmm_custom`
- Included custom-family constructors for Gaussian identity, binomial probit, binomial complementary-log-log, and quasi-Poisson log

### Adaptive quadrature

- Scalar grouped random coefficient through `fit_glmm_aghq`
- Correlated multidimensional random coefficients through `fit_glmm_aghq_multidimensional`
- Binomial, Poisson, and fixed-size negative-binomial responses
- Tensor-product adaptive Gauss-Hermite rules with a configurable node limit

### Nonlinear mixed models

- Gaussian nonlinear mixed models through `fit_nlmm`
- User-supplied nonlinear mean callback
- One grouped vector of random coefficients
- Laplace-integrated random effects and BOBYQA outer optimization
- Prediction and simulation

### Inference and grouped utilities

- Wald confidence intervals
- Profile-likelihood confidence intervals for fixed effects
- Parametric bootstrap for LMMs and built-in GLMMs
- Percentile bootstrap intervals
- Group-deletion influence diagnostics and Cook distances
- Likelihood-ratio tests
- `lmList`-style grouped weighted regressions and prediction
- Covariance conversion, random-effect PCA, singularity checks, residual diagnostics, and simulation
- Standard-normal Gauss-Hermite rules through order 200

## Basic use

```fortran
use lme4

type(random_term_t) :: terms(1)
type(lmm_result_t) :: fit

terms(1)%z = random_design
terms(1)%group = group_index
terms(1)%n_levels = n_levels
terms(1)%name = 'subject'

call fit_lmm_pls(y, x, terms, fit, reml=.true.)
if (.not. fit%converged) error stop fit%message
```

A correlated random-intercept/random-slope binomial model with three-node tensor AGHQ is fitted by:

```fortran
call fit_glmm_aghq_multidimensional(y, x, term, family_binomial, fit, order=3)
```

A nonlinear mixed model supplies a callback with signature:

```fortran
function mean_function(covariates, beta, random_effect) result(mean)
   use lme4, only : dp
   real(dp), intent(in) :: covariates(:), beta(:), random_effect(:)
   real(dp) :: mean
end function mean_function
```

and is fitted by:

```fortran
call fit_nlmm(y, covariates, group, n_levels, n_random, mean_function, start_beta, fit)
```

## Build with FPM

```text
fpm build
fpm test
fpm run --example lme4_example
fpm run --example glmm_extensions_example
fpm run --example advanced_algorithms_example
```

The bundled dependency is registered as `minqa`, matching `dependencies/minqa/fpm.toml`.

## Scope and scalability

The new Woodbury/PLS LMM solver avoids forming an `n by n` marginal covariance matrix and is preferable when the number of random coefficients is much smaller than the number of observations. It still stores the expanded random-effect design and its factorizations densely. R `lme4` uses sparse Eigen/CHOLMOD machinery and remains substantially more scalable for large crossed or highly sparse models.

Multidimensional AGHQ currently supports one grouped vector-valued random-effect term. A tensor rule has `order**q` nodes per group, so high dimension is intentionally guarded by `max_nodes`.

See `docs/API_MAP.md` and `docs/PORTING_NOTES.md` for exact coverage and remaining limitations.

## License

The original package is GPL version 2 or later. This port is distributed under the same terms. The bundled `minqa` translation is GPL-2.0-only. Original source files are retained under `original_source/` for attribution and traceability.
