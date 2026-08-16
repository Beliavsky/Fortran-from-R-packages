# Changelog

## 0.1.0

- Initial modern Fortran/FPM port of DiscreteDists 1.1.2.
- Ported all 14 d/p/q/r distribution families.
- Added Fortran `discrete_family_t` analogue of computational GAMLSS family
  behavior.
- Ported all exported estimators/log-likelihood helpers and numerical helpers.
- Reused COMPoissonReg-fortran as the CMP runtime dependency.
- Corrected the DLD CDF support shift and standardised quantile `log_p`
  semantics.
- Omitted plotting-only `plot_discrete_cdf` and R/S3/formula infrastructure.
