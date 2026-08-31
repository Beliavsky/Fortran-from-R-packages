# jomo — modern Fortran computational translation

This directory is a modern free-form Fortran translation of the reusable
computational core of the R package **jomo 2.7-6**, *Multilevel Joint Modelling
Multiple Imputation*, by Matteo Quartagno and James Carpenter.

The translation targets FPM and gfortran and contains no R runtime, plotting,
formula parser, S3 interface, BLAS/LAPACK system link, or vendored translated
R-package dependency.  It uses a deterministic package-local random-number
state so tests are reproducible.

## Implemented numerical areas

- Single-level continuous, categorical latent-normal, and mixed joint-model MCMC.
- Homogeneous multilevel joint models with arbitrary random-effect design
  columns and a full joint random-effect covariance across outcomes/slopes.
- Cluster-specific level-1 covariance models, including the hierarchical
  inverse-Wishart scale draw and jomo's sampled inverse-Wishart degrees
  parameter `a`.
- Two-level models with continuous/categorical/mixed variables at both levels,
  including the joint covariance of random effects and level-2 residuals.
- Direct two-level heterogeneous (`jomo2hr`) models combining that joint
  structure with cluster-specific level-1 covariance and sampled hierarchy `a`.
- Conditional multivariate-normal missing-data draws.
- Constrained latent-normal covariance Metropolis updates for categorical data.
- Matrix-normal, multivariate-normal, Wishart and inverse-Wishart draws.
- Substantive-model-compatible (SMC) imputation for Gaussian, binary probit,
  ordinal probit, and ordered-risk-set Cox models.  The SMC driver reproduces
  the upstream proposal logic for substantive predictors, exact conditional
  Gaussian draws for auxiliary variables, and design rebuilding after every
  proposal for powers, categorical dummies, and interactions.
- Multilevel SMC substantive models with random intercepts/slopes and a full
  random-effect covariance, plus level-2 predictor proposals and optional
  cluster-specific level-1 covariance matrices for the heterogeneous branches.
- Substantive latent-response, threshold, Gaussian variance, random-effect
  covariance, and upstream-style Cox coordinate-Newton parameter updates.

See [`API_COVERAGE.md`](API_COVERAGE.md) for exact scope and remaining parity
work.

## Build and test

From this directory:

```text
fpm build
fpm test
fpm run --example jomo1_continuous_example
fpm run --example jomo_smc_example
```

No separately installed BLAS or LAPACK library is required by this package.

## Minimal use

```fortran
use jomo, only : dp, i8, rng_state, rng_seed, jomo1_result, jomo1con_mcmc
```

Arrays and masks are supplied directly.  Categorical outcomes use integer
levels `1..K`; missingness is represented by explicit logical masks rather than
R `NA` payloads.  Cluster labels are one-based and contiguous.  SMC formula
terms are represented by `smc_design_spec`: a term contains one or more factors,
so polynomial powers and interactions are represented without an R formula
parser.  Categorical factors expand to `K-1` indicators with level `K` as the
latent-normal reference category.

## Licensing and provenance

The R package is GPL-2.  This translation is distributed under
GPL-2.0-only.  The upstream native sources also contain helper routines marked
GNU LGPL and attributed to Barry Brown, James Lovato, Guannan Zhang, and John
Burkardt; those notices are retained in `NOTICE.md` and
`UPSTREAM_PROVENANCE.md`, and an LGPL-2.1 license text is included for
conservative notice preservation.  The Fortran random/distribution/linear
algebra implementations here are newly written and do not vendor the upstream
C helper source.
