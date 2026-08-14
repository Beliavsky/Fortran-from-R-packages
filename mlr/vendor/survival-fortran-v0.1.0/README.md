# survival-fortran

Modern Fortran/FPM translation of the core numerical algorithms from the R
package **survival 3.8-9** (Terry M. Therneau et al.).

The public Fortran API operates directly on numeric arrays rather than R
formulas, model frames, factors, S3/S4 objects, environments, or graphics.

## Implemented computational families

- Kaplan-Meier survival and Nelson-Aalen cumulative hazard for right-censored data.
- Start/stop counting-process Kaplan-Meier/Nelson-Aalen risk sets.
- Aalen-Johansen multi-state transition-probability recursion.
- Cox proportional-hazards fitting with Breslow or Efron ties.
- Cox fitting for start/stop counting-process data.
- Cox baseline cumulative hazard, martingale residuals, and Schoenfeld residuals.
- Parametric AFT (`survreg`-style) fitting for extreme-value, Weibull,
  exponential, Rayleigh, Gaussian/lognormal, and logistic/loglogistic families.
- Fleming-Harrington/log-rank `survdiff` statistics with optional strata.
- Harrell-style right-censored concordance counts and C-index.
- Low-level Fine-Gray interval expansion corresponding to `src/finegray.c`.
- `survSplit`-style interval splitting.
- Jackknife survival pseudo-values.
- `pspline` basis and second-difference penalty construction using the supplied
  `splines-fortran` translation.

See `TRANSLATION_COVERAGE.md` for exact scope and differences.

## Build

```text
fpm build
fpm test
```

The root package has a local FPM dependency on
`vendor/splines-fortran-v0.1.0`.

## Licensing

The translated `survival`-derived sources retain the upstream declaration
**LGPL (>= 2)** (`LGPL-2.0-or-later`). The supplied `splines-fortran`
dependency is **GPL-2.0-or-later**. Because the default FPM target links the
GPL spline dependency, the combined executable/library distribution is marked
GPL-2.0-or-later while the individual survival-derived source files retain
their LGPL notices. See `UPSTREAM_PROVENANCE.md`.
