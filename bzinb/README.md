# bzinb-fortran

Modern free-format Fortran/FPM translation of the computational code in R package
`bzinb` 1.0.8.

The library implements bivariate Poisson (BP), two bivariate zero-inflated
Poisson models (BZIP-A and BZIP-B), bivariate negative binomial (BNB), and
bivariate zero-inflated negative binomial (BZINB) models.

## Build

```text
fpm build
fpm test
fpm run --example demo_bzinb
```

Direct GNU Fortran validation scripts are also supplied:

```text
scripts/validate.sh
scripts\validate.bat
```

## Main API

Use the facade module:

```fortran
use bzinb
```

Probability functions:

- `bp_pmf`, `bp_logpmf`
- `bzip_a_pmf`, `bzip_a_logpmf`
- `bzip_b_pmf`, `bzip_b_logpmf`
- `bnb_pmf`, `bnb_logpmf`
- `bzinb_pmf`, `bzinb_logpmf`

Simulation:

- `rbp_sample`
- `rbzip_a_sample`
- `rbzip_b_sample`
- `rbnb_sample`
- `rbzinb_sample`

Fitting and likelihood:

- `fit_bp`, `loglik_bp`
- `fit_bzip_a`, `loglik_bzip_a`
- `fit_bzip_b`, `loglik_bzip_b`
- `fit_bnb`, `fit_bnb_em`, `loglik_bnb`
- `fit_bzinb`, `fit_bzinb_em`, `loglik_bzinb`
- `bzinb_standard_errors`

`fit_bnb` and `fit_bzinb` now use the package's specialized EM machinery by
default. The v0.1 transformed-likelihood alternatives remain available as
`fit_bnb_direct` and `fit_bzinb_direct`.

The EM fit result types expose:

- natural-scale parameter estimates;
- source score-outer-product information matrices;
- covariance matrices and standard errors;
- underlying `rho` and `logit(rho)` with delta-method SEs;
- best log likelihood and historical-best iteration;
- likelihood trajectory and convergence status.

Other translated computations:

- `idigamma` and `inverse_digamma`
- `true_correlation`
- `nondropout_weight`
- `weighted_pearson_correlation`
- `pairwise_bzinb`
- `pairwise_bzinb_full`

`pairwise_bzinb_full` provides feature-pair indices, rho/SE, nonzero
proportions, optional full parameter/SE records, log likelihoods, iteration
counts, and optional random pair subsampling.

## Specialized EM implementation

The main v0.2.0 change is a direct Fortran translation of the numerical engine in
upstream `src/expt.cpp`, `src/opt.cpp`, and `src/em.cpp`. It includes latent
expectation sums, the source rescaling rules for extreme intermediate values,
inverse-digamma shape updates, the coupled scale update, zero-inflation class
updates, historical likelihood maximum selection, and the upstream analytic
score outer-product information matrix.

The information matrix uses the eight free coordinates
`(a0,a1,a2,b1,b2,p1,p2,p3)`. `p4` is constrained to one minus the other three
mixture probabilities. For convenience, the Fortran result expands the inverse
information matrix to a 9 x 9 covariance matrix including `p4`.

See `docs/TRANSLATION_STATUS.md` for the small source quirks intentionally
preserved and one pathological C++ undefined-memory branch that is safely
regularized in Fortran.

## Numerical dependencies

The compiled package is dependency-free. Positive-argument digamma/trigamma,
inverse digamma, Gamma/Poisson simulation, optimization, and matrix inversion
are provided natively.

## License

The upstream package declares `GPL-2`. This translation is distributed under
GPL-2.0-only. The complete upstream source is preserved under `upstream/`.
