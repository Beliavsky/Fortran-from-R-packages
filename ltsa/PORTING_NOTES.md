# Porting notes

## Data and object model

R vectors and matrices map to allocatable `real(dp)` arrays. R lists map to
`dl_ar_result`, `forecast_result`, `exact_likelihood_result`, and
`innovation_variance_result`. Errors are returned in `type(ltsa_error)` instead
of calling `stop()`.

R random-generator callbacks are represented by either the library's seeded
normal generator or an explicitly supplied innovation array. This keeps the
numerical routines deterministic and avoids procedure-pointer state in the
main API.

## Algorithms retained

- Durbin-Levinson coefficient, PACF, residual, likelihood, and simulation
  recursions are direct free-form Fortran translations.
- The Trench inverse uses the upstream Durbin-Levinson/upper-wedge construction.
  A checked Cholesky inverse is used only if a floating-point residual test
  detects loss of accuracy.
- The bordering formula in `ToeplitzInverseUpdate` is translated directly.
- `tacvfARMA` follows the upstream McLeod linear-system construction and its
  partial-autorrelation stationarity test.
- Davies-Harte simulation uses an in-package radix-2 complex FFT.

## Deliberate corrections and compatibility controls

1. The pure-R branch of upstream `DLSimulate` computes the first innovation but
   never assigns `z[1]`. The compiled C path does assign it. The Fortran port
   follows the compiled C behavior.
2. Upstream `DHSimulate` creates its two real Fourier endpoints as
   `2 + sqrt(2) * rnorm(2)`, introducing a nonzero mean. The Fortran default
   uses the standard zero-mean Davies-Harte construction. Set
   `source_compatible=.true.` to reproduce the upstream endpoint expression.
3. Upstream `exactLoglikelihood(..., innovationVarianceQ=FALSE)` returns `NA`
   despite documentation describing that mode. The Fortran routine implements
   the full known-covariance Gaussian likelihood instead.
4. The white-noise branch of upstream `tacvfARMA` refers to `maxLagp1` before it
   is assigned. The Fortran white-noise result is correctly returned as
   `(sigma2, 0, ..., 0)`.

## Numerical equivalents

The upstream AR innovation-variance estimator calls `stats::ar(...,
method="burg", aic=TRUE)`. To remain self-contained, the Fortran implementation
uses the biased sample autocovariance, Durbin-Levinson/Yule-Walker fits, and the
same AIC-style order tradeoff. The Kolmogoroff method uses a raw periodogram and
an optional centered moving-average span rather than reproducing every
`spec.pgram` taper and kernel option.

`SimGLP` uses direct convolution. This is algebraically identical to the
upstream FFT convolution and is preferable for modest filter lengths.
