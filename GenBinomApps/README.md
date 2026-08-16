# GenBinomApps-fortran

Modern Fortran/FPM translation of `GenBinomApps` 1.2.1.

Implemented:

- generalized/Poisson-binomial PMF, CDF, quantile, and RNG;
- ordinary one-sided and two-sided Clopper-Pearson confidence intervals;
- countermeasure-adjusted Clopper-Pearson confidence intervals;
- ordinary required sample-size calculation;
- countermeasure-adjusted required sample-size calculation;
- standalone regularized incomplete-beta and inverse-beta numerics.

No external runtime dependencies are required.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The original R package is retained under `upstream/`.
