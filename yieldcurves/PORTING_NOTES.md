# Porting notes

## Data structures

R lists, data frames, and S3 classes are represented by derived types:

- `curve_t`
- `series_t`
- `duration_result_t`
- `bond_duration_result_t`
- `zspread_result_t`
- `carry_result_t`
- `slope_result_t`
- `factor_result_t`
- `pca_result_t`

Each type carries an `ok` flag and diagnostic message where appropriate.

## Curve fitting

For a fixed decay parameter, Nelson-Siegel and Svensson are linear regression
models. The Fortran implementation profiles out the beta coefficients rather
than optimizing all parameters simultaneously. This is numerically more
stable and minimizes the same weighted sum of squared residuals.

Nelson-Siegel uses bounded scalar optimization on `0.01 <= tau <= 30`.
Svensson uses the upstream multi-start decay grid followed by bounded
Nelder-Mead optimization of `tau1` and `tau2`.

## Splines

The natural spline uses zero second derivatives at both endpoints. The FMM
implementation uses not-a-knot endpoint conditions, matching the numerical
behavior of R's `splinefun(method="fmm")` for ordinary data sets.

## Rate conversions

Coupon maturities must align with the selected annual or semi-annual coupon
frequency. The R package warns for some non-aligned annual maturities; the
Fortran API returns a failed status because a principal cash-flow date would
otherwise be ambiguous.

## Z-spreads and key-rate durations

To preserve upstream behavior, fitted curves are sampled at their original
maturities and linearly interpolated to bond cash-flow dates. Z-spreads use
discrete annual discounting `(1 + zero_rate + spread)^(-time)`, as in the R
implementation. Key-rate durations use the same one-sided triangular bump.

## PCA

R's `prcomp` uses an SVD. The Fortran implementation diagonalizes the centered
sample covariance or correlation matrix. The resulting principal subspace,
scores, component standard deviations, and variance shares are equivalent;
component signs are normalized deterministically because eigenvector signs
are otherwise arbitrary.

## Missing values

The upstream functions reject missing values. Fortran inputs must be finite.
Unavailable standard slope tenors are represented by IEEE quiet NaNs.
