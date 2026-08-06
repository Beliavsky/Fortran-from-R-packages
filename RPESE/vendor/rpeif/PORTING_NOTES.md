# Porting notes

## Scope

The translation covers every exported computational influence-function family,
nuisance-parameter construction, robust cleaning, shape-grid generation, and
optional AR prewhitening. R plotting and time-series class infrastructure are
not translated.

## Numerical conventions

### Quantiles

The empirical quantile implementation uses R's default type-7 interpolation.

### Kernel density at a quantile

RPEIF calls `density()` and interpolates the result. This port uses the same
`bw.nrd0` bandwidth rule but evaluates the Gaussian kernel estimate directly at
the requested point. This avoids an FFT grid and interpolation. Results should
be close but are not expected to be bit-for-bit identical to `stats::density`.

### Prewhitening

The R package uses `stats::arima` for an AR(p) model with a mean. The Fortran
port uses conditional least squares with an intercept. The first `p` residuals
are set to zero so the returned vector has the original length. This is a
self-contained computational equivalent, not an exact maximum-likelihood ARIMA
replication.

### Robust mean and cleaning

The psi, psi-derivative, weight, scale, and efficiency-tuning functions are
reused from the completed RobStatTM translation. Location is found by iterated
reweighted averaging starting at the median with MAD scale. Cleaning clips at
the upstream constants for efficiencies 0.95, 0.99, and 0.999. For other valid
efficiencies, a normal-quantile cutoff is used rather than leaving the cutoff
undefined.

## Source compatibility

`rpeif_options%source_compatibility` defaults to `.true.`.

### Omega UPM sign

The upstream `UPM` helper computes `(const - return)^order` for observations
above the threshold. At order 1 this is negative, while the nuisance-parameter
formula and standard UPM definition are positive. Source-compatible mode keeps
the negative empirical value. Corrected mode uses `(return - const)^order`.

### Empirical VaR-ratio denominator

The upstream empirical VaR-ratio routine assigns the empirical quantile to a
variable named `fq.alpha` and then uses it where the density at the quantile is
required. Source-compatible mode keeps this behavior. Corrected mode uses the
Gaussian KDE density, consistent with the VaR influence function and nuisance
formula.

### Sharpe-ratio risk-free rate

The upstream routine computes the Sharpe ratio with `mu-rf` but uses `mu` in
the scale contribution to the influence function. Source-compatible mode keeps
that formula. Corrected mode uses `mu-rf` in that term.

### DSR nuisance shape

The theoretical nuisance path keeps the upstream stored `dsr` value in
source-compatible mode. Corrected mode recomputes it using the supplied
risk-free rate and the data-path definition.

## Error handling

R exceptions are replaced by integer status codes and explanatory messages.
Numerically undefined influence functions return zero-filled values and
`rpeif_numerical_failure` rather than NaN or infinity.

## Tail probabilities and nuisance parameters

The upstream default nuisance object uses `alpha=0.1`, while the VaR, ES, and
VaR-ratio wrappers default to 0.05. The Fortran port preserves this separation:
formula probabilities come from `rpeif_options`, while quantiles and tail
expectations in theoretical evaluation come from the supplied nuisance object.
For internally consistent theoretical shapes, construct the nuisance object
with the same alpha and beta used in the options.
