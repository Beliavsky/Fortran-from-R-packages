# ordinal

Modern free-form Fortran translation of the computational core of the R package
**ordinal** 2026.7-26 by Rune Haubo Bojesen Christensen.

The library implements cumulative-link models for ordered responses, including
nominal and scale effects, plus Gaussian cumulative-link mixed models. It is a
numerical API: callers provide integer category codes and already-constructed
numeric design matrices. R formula parsing, S3 methods, plotting, and `Matrix`
object machinery are deliberately outside the translation scope.

## Build

```text
fpm build
fpm test
fpm run --example clm_example
fpm run --example clmm_general_example
```

No system BLAS/LAPACK installation is required. The implementation is
self-contained and does not vendor or copy `numDeriv`, `nlme`, BLAS, LAPACK, or
other translated packages.

## Main modules

- `ordinal`: public facade, including `dp`.
- `ordinal_links`: link CDFs, densities, and density derivatives.
- `ordinal_thresholds`: flexible/symmetric/symmetric2/equidistant cut points.
- `ordinal_clm`: CLM likelihood, analytic derivatives, fitting, nominal/scale
  effects, and prediction.
- `ordinal_clmm`: scalar random-intercept CLMM with Laplace, adaptive GHQ, and
  non-adaptive GHQ paths following `nAGQ` semantics.
- `ordinal_clmm_general`: dense Laplace CLMM for multiple terms, random slopes,
  crossed/nested groups, and correlated random effects.
- `ordinal_profile`: CLM profile likelihood and profile/Wald intervals.
- `ordinal_rank`: rank detection and rank-deficient-column reduction.
- `ordinal_quadrature`: generated Gauss-Hermite quadrature rules.
- `ordinal_numerics`: deterministic numerical optimization, matrix helpers,
  finite differences, and Hessian diagnostics.

`X` matrices passed to CLM/CLMM initializers contain location predictors
**without an intercept**. Thresholds play the intercept role in cumulative-link
models. Nominal and scale design matrices likewise omit implicit intercepts
unless the caller intentionally wants one as a numeric column.

## Mixed-model scope

`init_clmm_problem` is the compact one-random-intercept API and supports the
upstream-style `nAGQ` choices. For random slopes or multiple random-effect terms,
use `init_clmm_laplace_problem`; columns of `re_z` are concatenated by term and
`term_q` gives the number of random coefficients in each term.

See `API_COVERAGE.md` for detailed parity notes and `NOTICE` for upstream
provenance.
