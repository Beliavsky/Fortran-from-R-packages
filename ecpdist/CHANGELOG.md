# Changelog

## 0.1.0

- Initial modern Fortran/FPM computational translation of `ecpdist` 0.2.1.
- Added ECP density, CDF, survival, quantile, hazard, cumulative hazard, and RNG.
- Added raw moments, conditional moments, mean residual life, Bowley skewness,
  and Moors kurtosis.
- Added stable log-survival probability evaluation for positive/negative `phi`.
- Added adaptive 15-point Gauss-Kronrod quantile integration.
- Corrected upstream `qecp` `lower_tail` and `log_p` option behavior.
- Omitted plotting-only `ecp_plot`.
