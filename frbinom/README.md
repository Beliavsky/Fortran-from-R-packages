# frbinom-fortran

Modern Fortran/FPM translation of the R package `frbinom` 1.0.0.

The package implements two fractional-binomial distributions based on
generalized Bernoulli processes with dependence.

Implemented:

- PMF, CDF, quantile, and RNG for fractional binomial family I;
- PMF, CDF, quantile, and RNG for fractional binomial family II;
- `start` variants for both families;
- scalar and vector evaluation;
- lower- and upper-tail CDF/quantile interfaces;
- reusable PMF/CDF table builders;
- deterministic Fortran RNG seeding helper.

No external runtime dependencies are required.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The complete supplied R package is retained under `upstream/`.
