# Porting notes

## Data representation

An R formula such as

```text
y ~ x + (1 + x | subject)
```

is represented by `x(:,1:2)`, `term%z(:,1:2)`, integer group identifiers, and an explicit level count. Multiple independent grouping terms are supplied as an array of `random_term_t` values.

The `covariance_structure` member selects:

- `covariance_unstructured`;
- `covariance_diagonal`;
- `covariance_compound_symmetry`;
- `covariance_ar1`.

Structured forms reduce the number of optimized covariance parameters and are constructed so the resulting covariance matrix remains positive definite within the parameter bounds.

## LMM likelihood formulations

`fit_lmm` uses a relative random-effect covariance `G` and prior weights `w` to form

```text
V0 = diag(1/w) + Z G Z'
```

and profiles the residual variance. ML and REML criteria include Gaussian constants, `log|V0|`, and for REML `log|X' V0^-1 X|`.

`fit_lmm_pls` evaluates the same likelihood with the Woodbury identity. It factors

```text
C = Z' W Z + G^-1
```

and computes products with `V0^-1` without constructing `V0`. This is substantially more memory efficient when the total random-effect dimension is smaller than the number of observations, but `Z`, `C`, and related matrices remain dense.

## General and custom GLMMs

For each covariance parameter vector, PIRLS jointly solves for fixed effects and conditional random-effect modes. The Laplace objective combines the conditional negative twice-log-likelihood, the random-effect quadratic penalty, `log|G|`, and `log|Z'WZ + G^-1|`.

`family_spec_t` replaces R family closures with procedure pointers. A custom family supplies scalar procedures for:

1. link;
2. inverse link;
3. derivative of the mean with respect to the linear predictor;
4. variance;
5. per-observation log likelihood;
6. response validation.

The caller is responsible for numerical consistency among those procedures and for defining valid behavior near the response-domain boundary.

## Adaptive Gauss-Hermite quadrature

`fit_glmm_aghq` handles one grouped scalar random coefficient. `fit_glmm_aghq_multidimensional` generalizes the same construction to a vector random effect for one grouping factor:

1. find each group's posterior mode with Newton iteration;
2. factor the negative posterior Hessian at the mode;
3. transform a tensor product of standard-normal Gauss-Hermite nodes;
4. calculate the group marginal likelihood by log-sum-exp;
5. optimize fixed effects and covariance parameters jointly with BOBYQA.

The number of evaluations per group is `order**q`; `max_nodes` prevents accidental combinatorial explosion. Multiple independent or crossed random-effect terms are still fitted with the Laplace routine rather than joint AGHQ.

## Nonlinear mixed models

`fit_nlmm` takes a nonlinear mean callback `mean(covariates,beta,random_effect)`. It supports a Gaussian residual model and one grouped random-effect vector. Group modes and Hessians are calculated using finite differences, integrated with a Laplace approximation, and nested inside bounded BOBYQA optimization of fixed effects, covariance parameters, and residual scale.

This is the numerical core of a useful `nlmer`-style model but does not parse nonlinear formulas or symbolic derivative attributes.

## Inference

Wald intervals use the estimated fixed-effect covariance. Profile intervals constrain one fixed effect, remove its design column, adjust the response or offset, and refit all nuisance parameters until the likelihood-ratio target is reached.

Parametric bootstrap simulates from the fitted model and refits each sample. Failed refits are retained as non-finite entries and excluded by percentile interval calculation. Influence diagnostics delete one grouping level at a time and report DFBETA, Cook distance, and deleted log likelihood.

## Simulation

Gamma simulation uses Marsaglia-Tsang sampling. Negative-binomial simulation uses a Gamma-Poisson mixture. Inverse-Gaussian simulation uses the Michael-Schucany-Haas transformation. Gaussian nonlinear simulation draws grouped multivariate-normal random effects and Gaussian residuals.

## Remaining sparse limitation

The port has no dependency on SuiteSparse, Eigen, or CHOLMOD. Consequently it cannot reproduce sparse symbolic factorization, sparse rank updates/downdates, or the memory scaling of R `lme4`. Adding that capability would require an external sparse backend and a separate interoperability layer rather than a self-contained translation.
