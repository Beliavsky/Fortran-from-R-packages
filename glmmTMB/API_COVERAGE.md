# API and computational coverage

This document distinguishes translated numerical functionality from upstream
R/TMB runtime machinery that is intentionally outside this standalone Fortran
port.

## Response families

| Upstream family | Fortran conditional log likelihood | Notes |
|---|---:|---|
| Gaussian | yes | `mu`, data-scale `phi` standard deviation |
| skew-normal | yes | Upstream mean/SD/slant reparameterization |
| Poisson | yes | Stable log PMF |
| binomial | yes | Uses logit-scale robust binomial kernel |
| Gamma | yes | Shape=`phi`, scale=`mu/phi` |
| beta | yes | Ferrari-Cribari-Neto mean/precision form |
| ordered beta | yes | Boundary masses plus beta interior |
| beta-binomial | yes | Robust log-shape implementation |
| nbinom1 | yes | `Var=mu*(1+phi)` |
| nbinom2 | yes | `Var=mu*(1+mu/phi)` |
| nbinom12 | yes | Extra `psi` dispersion component |
| truncated nbinom1 | yes | Exact zero-probability correction |
| truncated nbinom2 | yes | Exact zero-probability correction |
| truncated Poisson | yes | Exact zero-probability correction |
| generalized Poisson | yes | Upstream theta/lambda mapping |
| truncated generalized Poisson | yes | Exact zero-probability correction |
| COM-Poisson | yes | Exact-mean numerical lambda solve |
| truncated COM-Poisson | yes | Exact zero-probability correction |
| Tweedie | yes | 1<p<2 compound-Poisson/gamma series |
| lognormal | yes | Parameterized by mean and SD on data scale |
| Student-t | yes | `df=exp(psi(1))` |
| Bell | yes | Lambert-W back transformation |

Zero inflation is implemented with the same log-space mixture calculation as the
upstream C++ objective.  Observation-level NaNs are skipped by the vector NLL
routine, corresponding to upstream missing-response behavior.

## Link functions

Forward/inverse link support includes log, logit, probit, inverse,
complementary-log-log, identity, and square-root links.  The upstream Bell
Lambert-W inverse link is included.  The upstream C++ source does not implement
the forward Lambert-W link, and neither does this translation.

Stable helpers corresponding to `logit_inverse_linkfun`,
`log_inverse_linkfun`, `log1m_inverse_linkfun`, and `calc_log_nzprob` are
translated.

## Random-effect covariance structures

`covariance_term_nll` implements:

- diagonal (`diag`);
- homogeneous diagonal (`homdiag`);
- unstructured (`us`);
- compound symmetry (`cs`) and homogeneous compound symmetry (`homcs`);
- Toeplitz (`toep`) and homogeneous Toeplitz (`homtoep`);
- AR(1) and heterogeneous AR(1);
- continuous-time Ornstein-Uhlenbeck;
- exponential spatial correlation;
- Gaussian spatial correlation;
- Matérn spatial correlation;
- reduced-rank factor loading (`rr`);
- proportional unstructured (`propto`);
- equal-to unstructured (`equalto`).

The negative log density and implied correlation/standard-deviation outputs are
available.  Reduced-rank factor loadings are also available.  Upstream random
simulation branches are excluded because they are wired to R/TMB RNG services.

## Distribution/prior helpers

Translated helpers include generalized Poisson, robust beta-binomial,
skew-normal, Cauchy, Bell/Lambert-W, COM-Poisson, Tweedie, multivariate log
gamma, LKJ, Wishart, and inverse-Wishart kernels.

Normal, gamma, Student-t, Cauchy, and LKJ priors used by the upstream objective
are implemented.  `beta_prior` is declared in the upstream enum but is not
implemented by the upstream C++ prior switch; this translation preserves that
status rather than inventing behavior.

The upstream gamma prior evaluates a gamma density at `exp(parameter)` without
adding a log-Jacobian term.  `scalar_prior_log_density` deliberately preserves
that exact behavior.

## Model composition

- `build_linear_predictor` computes dense `X*beta + Z*b + offset`.
- `observation_nll_vector` sums weighted conditional observation likelihoods.
- `random_terms_nll` combines split Gaussian random-effect terms.
- `glmmtmb_joint_nll` combines observation likelihood, Gaussian random-effect
  density, and a supplied prior NLL.

These routines expose the joint objective used before random effects are
integrated out.

## Intentionally not translated

The following are not represented as standalone Fortran APIs:

- TMB automatic-differentiation tape construction and reverse-mode gradients;
- Laplace integration of random effects and the optimizer that minimizes the
  resulting marginal objective;
- TMB/R sparse-matrix wrappers and RcppEigen object plumbing;
- R formula parsing, contrasts, model frames, factor handling, S3/S4 methods,
  prediction object reconstruction, printing, plotting, diagnostics, and
  package-registration code;
- R parallel-cluster/OpenMP control interfaces;
- R RNG-backed simulation branches;
- denominator-degree-of-freedom/R model-object orchestration through packages
  such as `pbkrtest`;
- serialization of fitted R model objects.

Consequently, this is substantial computational parity for the portable joint
likelihood and covariance/family kernels, but not end-to-end parity with the R
`glmmTMB()` user interface or its AD/Laplace fitting engine.

## Known numerical differences

1. Upstream COM-Poisson normalization/mean inversion uses TMB's specialized
   atomic implementation.  The Fortran version uses log-scale finite summation
   and bisection and intentionally rejects extremely large support requirements.
2. Upstream Matérn correlation uses TMB's arbitrary-order Bessel K support.  The
   Fortran implementation uses midpoint quadrature of the integral
   representation of K_nu.
3. Upstream Tweedie uses TMB's `dtweedie`; the Fortran translation evaluates the
   mathematically equivalent 1<p<2 compound-Poisson/gamma series directly.
4. `family_variance` reproduces the upstream R skew-normal variance expression,
   even though the C++ skew-normal likelihood parameterizes `phi` as a
   data-scale standard deviation.  This upstream inconsistency is preserved and
   documented rather than silently changed.
