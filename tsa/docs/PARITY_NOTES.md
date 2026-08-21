# Exact-parity notes for v0.5.0

## TSA `arimax` transfer blocks

TSA increments each requested transfer MA order by one before fitting. A user
order `(p,q)` therefore estimates `p` recursive denominator coefficients and
`q+1` convolution numerator coefficients for lags `0:q`. The Fortran
`transfer_spec` keeps the user-facing order and applies the increment internally.

Transfer contributions are recursive-filtered by their denominator and then
one-sided-convolved by their numerator before subtraction from the response.
Innovational-outlier pulse columns are rebuilt from the current trial ARIMA
parameters. Guarded filtering maps overflow-prone, unstable optimizer trial
regions to a large finite objective so strict IEEE traps do not abort fitting.

## ARIMA ML state space

The ML path now follows the state construction used by R `stats::arima`: the
original series is not pre-differenced, differencing coefficients enter the
observation/state equations, integrated components receive diffuse `kappa`
variance, and missing observations perform prediction-only updates. Diffuse
innovations with very large gain are omitted from the concentrated likelihood.
Returned ML residuals are standardized innovations.

The stationary ARMA block uses a direct translation of the default
Gardner/AS154 `C_getQ0` initialization. A discrete Lyapunov solver is retained
only as a large-order/allocation fallback.

## Optimization

TSA calls R `optim` with BFGS. v0.3.0 adds the corresponding variable-metric
BFGS path with optimization in scaled coordinates and central finite-difference
derivatives using the R default scaled perturbation. ARMA parameter scales are
one. Regression and transfer scales follow TSA's `10 * standard error` rules.
v0.4.0 additionally rotates multiple all-free regression columns by direct
right singular vectors, maps coefficients/covariances back to the original
basis, and implements R's numerical `optimHess` derivative-of-gradient scheme.
The stationarity covariance transform uses TSA's analytic Levinson derivative.
The remaining low-level numerical difference is that matrix factorization is a
self-contained pivoted LU rather than the platform LAPACK used by R.

## Regression SVD conditioning

When more than one regression coefficient is free, R/TSA rotate the complete
regressor matrix with the right singular vectors before optimization. v0.4.0
uses a direct one-sided Jacobi SVD and performs the fit in that basis. Reported
coefficients and covariance are transformed back with the same orthogonal
matrix. If any regression coefficient is fixed, no rotation is performed.

For TSA's specialized IO path, the R source rotates the matrix but retains the
old column names and then detects IO columns by name. That can filter a rotated
mixture rather than the original pulse. The Fortran port applies dynamic IO
filtering to the intended original IO effects and rotates the resulting design;
this is an intentional model-preserving correction, not literal reproduction
of that naming quirk.

## Spectra delegated by TSA

TSA's `spec` wrapper delegates to R rather than carrying its own spectral
implementation, so the current R `stats` algorithms are the parity reference.
`spec_pgram` implements the univariate and matrix periodogram paths, compact or
full symmetric smoothing kernels, named Daniell/modified-Daniell/Fejer/Dirichlet
kernels, scalar or per-series tapers, cross-spectra, squared coherence, and
phase. `spec_ar` now supports Yule-Walker, Burg/Burg2, OLS, and ML fitting with
fixed-order or AIC order selection. The Burg update is translated directly from
R's current C implementation, including both innovation-variance sequences. The ML method uses the translated modern
ARIMA state-space engine rather than R's legacy `arima0`, so tiny numerical
differences can remain there. Plotting remains outside scope.

## Bootstrap

The original `arima.boot` helper uses the fitted MA vector itself as the
convolution filter rather than prepending a unit coefficient. The Fortran port
preserves that source behavior. For conditional bootstrap, it preserves the
first `p+d` observations and generates the remainder.

## TAR matrices and missing values

TSA stacks complete rows from all input series and uses a baseline coefficient
block plus zero-filled series-specific deltas. The Fortran multiseries fitter
uses that layout. IEEE NaN plays the role of R `NA`; incomplete lag windows are
omitted. The documented centering/standardization operations are implemented
rather than reproducing the two apparent source bugs described in
`TRANSLATION_STATUS.md`.
