# mitml

Modern free-form Fortran translation of the reusable computational routines in
R package **mitml 0.4-5**, "Tools for Multiple Imputation in Multilevel
Modeling".

Upstream mitml was authored by Simon Grund, Alexander Robitzsch, and Oliver
Luedtke and is distributed under `GPL (>= 2)`. This translation retains that
license as `GPL-2.0-or-later`. See `NOTICE.md`, `PROVENANCE.md`, and `LICENSE`.

## Scope

The Fortran API concentrates on numerical work that is useful independently of
R's model objects and data-frame/S3 infrastructure:

- Rubin pooling of multiple-imputation estimates and covariance matrices;
- Barnard-Rubin finite-complete-data degrees of freedom, RIV, FMI, t tests, and
  confidence intervals;
- D1 and D2 Wald-type MI tests;
- D3 and D4 likelihood-ratio MI tests;
- linear-constraint testing and a transformed-estimand entry point for general
  delta-method constraints;
- NaN-aware cluster means, nested group/cluster means, and leave-one-out means;
- multilevel RB1, RB2, SB, MVP R-squared measures and ICC;
- Gelman-Rubin Rhat, SD-proportion/effective-sample-size diagnostics, reduced
  autocorrelation, and centered moving averages;
- Gaussian `lm` and `lmm` log-likelihood kernels used by mitml model-comparison
  calculations.

The upstream `panImpute()` and `jomoImpute()` functions are predominantly R
formula/data-management front ends. Their underlying numerical engines are
available as separate top-level translations `pan` and `jomo`. They are not
linked simultaneously into this package: the translated `pan` package is
GPL-3-only and the translated `jomo` package is GPL-2-only, so a single static
FPM dependency graph containing both would create an avoidable GPL-version
compatibility problem.

## Dependencies

Place this directory at the repository root alongside:

- `rfortran-core` -- supplies the shared `dp` kind and R-compatible t/F
  distribution functions;
- `rfortran-linalg` -- supplies explicit-interface dense linear solves and
  log determinants for Gaussian mixed-model likelihoods.

The FPM manifest uses sibling path dependencies and does not vendor their
source.

## Build and test

```text
fpm build
fpm test
```

An example is available through FPM's example target.

This source does not require a separately installed BLAS or LAPACK library.
`rfortran-linalg` uses the repository's pure-Fortran linear-algebra stack.

## Public module

Use `module mitml`. The package re-exports `dp` from `r_kinds` and the status
constants/result types used by the numerical API.

Typical pooling call:

```fortran
use mitml, only : dp, pool_estimates, pooled_estimates

type(pooled_estimates) :: result
real(dp) :: qhat(2, 5)
real(dp) :: uhat(2, 2, 5)

! Fill qhat and uhat, then:
call pool_estimates(qhat, result, uhat, df_complete=120.0_dp)
```

See `API_COVERAGE.md` for the R-to-Fortran mapping and intentional exclusions.
