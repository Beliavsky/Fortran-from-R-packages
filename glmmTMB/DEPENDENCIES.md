# Dependency review

Before implementing dependencies, the public
`Beliavsky/Fortran-from-R-packages` repository was checked for reusable
translations and shared numerical modules.

Relevant translated packages/shared layers observed include `lme4`, `Matrix`,
`numDeriv`, `nlme`, `mgcv`, `sandwich`, `tweedie`, and the repository's focused
`rfortran-*` numerical layers.  Those packages are dependencies of the R-facing
upstream package, but their APIs are not required by the portable numerical
kernel implemented here once formula/model-object/data-frame orchestration is
excluded.

The only FPM dependency used by this package is:

```toml
TMB = { path = "../TMB" }
```

This directly reuses the previously translated TMB density and covariance
primitives.  No dependency source is copied into `glmmTMB/`.

The translated `tweedie` package was deliberately not introduced as a new
runtime dependency.  glmmTMB needs a narrow Tweedie density with its own
1<p<2 parameterization, and the local compound-Poisson/gamma implementation
avoids unnecessary compatibility/system-library dependencies.
