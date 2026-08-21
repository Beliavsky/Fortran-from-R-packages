# Translation status

## Implemented through v0.5.0

The package translates the principal numerical/statistical algorithms from TSA
1.3.1. Versions 0.2.0 through 0.5.0 concentrate on source-level parity for the
areas TSA delegates to R or implements through internal helper closures.

### Core TSA computations

- descriptive skewness and excess kurtosis using TSA's population-moment
  definitions;
- ACF, autocovariance, PACF, CCF, periodogram, harmonics, seasonal indices, and
  exact runs calculations;
- OLS AR fitting/order selection, ARMA spectra, prewhitening, EACF, Box-Cox AR
  profiling, and ARMA subset selection;
- Ljung-Box/Box-Pierce, McLeod-Li, Keenan, and Tsay tests;
- additive and innovational outlier statistics;
- QAR, GARCH, and TAR simulation;
- univariate and multiseries SETAR fitting, prediction, skeleton/cycle
  calculations, and threshold likelihood-ratio inference;
- generalized GARCH Box diagnostics through the supplied tseries port.

### v0.2.0 parity work

- rational `xtransf` transfer blocks and dynamic `io` interventions;
- `fixed`/`init` NaN masks and separate CSS, CSS-ML, and ML modes;
- correct AR/SAR stationarity and MA/SMA invertibility handling;
- TSA bootstrap conditioning and source-level MA convolution behavior;
- multiseries TAR stacking, missing-window omission, and preprocessing.

### v0.3.0 parity work

1. **Diffuse integrated ARIMA ML.** ML now keeps the original regression-
   adjusted series and represents ordinary/seasonal differencing in the state
   vector. Integrated states receive `kappa` diffuse variance, missing
   observations receive prediction-only updates, high-gain diffuse innovations
   are excluded from the concentrated likelihood, and returned ML residuals are
   standardized innovations. CSS `n.cond` and ML `n.cond=0` follow R/TSA
   semantics.
2. **Gardner/AS154 stationary initialization.** The default stationary ARMA
   covariance is computed with a direct Fortran translation of R's `C_getQ0`
   packed AS154 algorithm. For unusually high state orders or allocation
   failure, an algebraically equivalent discrete Lyapunov solver is retained as
   a bounded-memory fallback.
3. **BFGS and TSA parameter scaling.** The bundled optimizer now includes the
   variable-metric BFGS algorithm used by R `optim(method="BFGS")`, with scaled
   coordinates and R-style finite-difference steps. ARIMAX uses BFGS by
   default. Regression coefficients receive TSA's `10 * OLS standard error`
   scale, and each transfer block receives ten times the simple-regression
   slope SE; ARMA parameters retain unit scale.
4. **`spec.pgram`.** The Fortran spectral module implements R's detrending,
   cosine taper, taper corrections, zero padding/2-3-5 fast length, modified
   Daniell kernels, circular smoothing, kernel df/bandwidth, and one-sided
   scaling. Matrix input additionally returns auto-spectra, squared coherence,
   and phase in R's pair ordering.
5. **`spec.ar` default path.** The univariate Yule-Walker path includes R's
   default order search/AIC convention, prediction-variance correction, and
   frequency-domain spectrum formula.

### v0.4.0 ARIMA/ARIMAX regression and covariance parity

1. **SVD-conditioned regressors.** Multiple regression columns are rotated by
   direct right singular vectors whenever all regression coefficients are free,
   matching the R/TSA conditioning rule. Optimization occurs in the rotated
   basis; reported coefficients and covariance matrices are mapped back to the
   original regressor basis. Fixed regression coefficients disable the rotation.
2. **R `optimHess` finite differences.** Covariance estimation now uses the
   derivative-of-scaled-gradient Hessian construction from R's `C_optimhess`
   rather than a generic second-difference Hessian.
3. **Analytic stationarity Jacobian.** The AR/SAR covariance transformation now
   uses the analytic derivative of TSA's Levinson/PACF recursion instead of a
   numerical Jacobian.
4. **Pivoted covariance solve.** A reusable pivoted-LU factor/solve path replaces
   repeated Gauss-Jordan inversion for the ARIMAX Hessian. It mirrors the
   numerical structure of a general LAPACK solve while keeping the FPM package
   self-contained.
5. **Regression initialization/scaling.** ML starting regressions now use
   complete cases. Plain ARIMA/ARIMAX regressors follow the differenced
   `stats::arima` scale calculation; the specialized TSA transfer/IO branch
   retains TSA's raw-series scale calculation.

### v0.5.0 spectral parity

1. **All univariate `spec.ar` methods.** `spec_ar` now exposes Yule-Walker,
   Burg/Burg2, OLS, and ML fitting. Fixed-order and AIC order-selection paths
   follow the corresponding R `ar.*` conventions, including a direct
   translation of the current C Burg recursion, both Burg innovation-variance
   choices, and the OLS embedded-regression criterion. The ML
   route reuses the translated diffuse/state-space ARIMA likelihood.
2. **Compact `tskernel` representation.** R-style compact symmetric kernel
   coefficients (lag zero followed by positive lags) can be supplied directly
   and are expanded internally. Named constructors are provided for Daniell,
   modified Daniell, Fejer, and Dirichlet kernels; Daniell variants accept
   vector bandwidths and compose them by convolution as R does.
3. **General spectral smoothing.** `spec_pgram` accepts either compact kernel
   coefficients or full symmetric weights in addition to the existing `spans`
   interface. Kernel df and bandwidth are computed from the actual smoothing
   weights.
4. **Per-series tapering.** Multivariate periodograms accept one taper fraction
   per series. Auto-spectra and df use their own taper-energy corrections,
   while cross-periodogram scaling preserves the coherence calculation.

## Deliberately omitted

These presentation/infrastructure layers remain outside scope:

- `lagplot` and all other plotting functions;
- S3 dispatch, formulas, `ts` attributes, printing, and summary formatting;
- R callback plumbing for arbitrary user-supplied RNG functions.

## Exact-parity work still remaining

1. **`ar.mle` legacy backend details.** R's `ar.mle` obtains candidate fits
   through its legacy `arima0` path. The Fortran `spec_ar(method="mle")` uses
   the package's newer translated diffuse/state-space ARIMA ML engine. The
   statistical model and order search are equivalent, but optimizer/startup
   details can produce small numerical differences.
2. **Linear-algebra backend bitwise differences.** SVD, pivoted LU, and related
   routines are self-contained Fortran implementations rather than R's
   platform BLAS/LAPACK. Ordinary floating-point agreement is expected, but
   nearly singular cases can differ in the last digits or singularity decision.
3. **Specialized IO + SVD source quirk.** TSA rotates the complete xreg matrix
   and then identifies IO columns by preserved names. The Fortran port instead
   filters the intended original IO effects before rotating the resulting
   design, preserving the intended model space rather than reproducing that
   name-dependent quirk.
4. **R runtime/object semantics.** IEEE-NaN handling is implemented where
   computationally material, but exact R recycling, attributes, warning text,
   `ts` metadata/indexing, optimizer call counts, arbitrary `tskernel` object
   attributes, and object-level NA behavior are outside the native Fortran API.

## Intentional corrections to upstream rough edges

Two `tar` arguments are documented by TSA but problematic in the R source:
`center` is calculated but not applied, and the logical argument `standard`
shadows the apparent standardization function. The Fortran port implements the
documented centering/standardization behavior rather than reproducing those
source bugs. See `docs/PARITY_NOTES.md`.
