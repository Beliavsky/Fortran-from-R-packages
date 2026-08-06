# Porting notes

## Numerical correspondence

The original fitting methods are preserved:

1. Construct the same candidate sequences with steps 0.5 (Nelson-Siegel), 1.0 and 1.5 (Svensson).
2. For each candidate target maturity, maximize the curvature loading over the same bounded interval.
3. Estimate the beta coefficients by linear least squares.
4. Select the candidate with the smallest sum of squared residuals.

The Fortran optimizer uses the R `optimize` default tolerance scale, `epsilon(1.0_dp)**0.25`.

## Intentional robustness changes

- Rank-deficient Svensson designs are rejected instead of allowing a partially identified regression to compete with zeroed coefficients.
- Maturities are required to be finite, positive, and strictly increasing.
- The small-argument factor formulas use Taylor expansions to avoid cancellation.
- IEEE NaN rate observations are omitted from regressions, corresponding to R's usual `na.omit` behavior.

## R infrastructure

`xts` timestamps and class restoration are not computational parts of the model. The Fortran matrix row order is the time order supplied by the caller. Plotting occurs only in R documentation examples and is not translated.
