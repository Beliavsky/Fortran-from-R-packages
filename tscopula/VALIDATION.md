# Validation

The package is tested in two configurations:

- Checked: GNU Fortran 2018, warnings as errors, bounds/shape/runtime checks,
  backtraces, and initialized invalid sentinel values.
- Optimized: GNU Fortran 2018 with `-O3` and warnings as errors.

The six test programs cover:

1. All margin families and V-transform inversion/probability identities.
2. Pair-copula CDFs, densities, rotations, Kendall conversions, and inverse
   h-functions.
3. ARMA/SARMA expansion, standardized variance, exact filtering, simulation,
   residuals, forecasts, and Kendall-PACF calculations.
4. D-vine simulation, likelihood, conditional density, Rosenblatt inversion,
   residuals, and independence reduction.
5. V-copula, W-copula, full-model, empirical-margin, forecasting, and AICc
   calculations.
6. ARMA, margin, D-vine, and full-model fitting plus compatibility APIs.

The example constructs an ARMA-generated Gaussian D-vine, applies a nonlinear
V-transform and skew-Student margin, simulates 500 observations, and reports
Kendall dependence, a one-step tail quantile, and copula likelihood.

Random sequences are reproducible within this library after `set_seed`, but
are not bit-for-bit identical to R's random-number streams.
