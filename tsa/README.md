# TSA-fortran

Modern free-format Fortran translation of the computational routines in the
R package **TSA 1.3.1** by Kung-Sik Chan and Brian Ripley. The project is laid
out for the Fortran Package Manager (FPM) and uses `dp = kind(1.0d0)`.

The original TSA package is preserved verbatim in `upstream/TSA-1.3.1/`.
The user-supplied Fortran translations of TSA's dependencies are preserved in
`vendor/`. The subset of the supplied `leaps` and `tseries` sources needed by
the translated library is copied under `src/dependencies/`, making the FPM
package self-contained.

## Build and test

```text
fpm build
fpm test
fpm run --example demo_tsa
```

The translation has also been compiled and tested directly with GNU Fortran
14 using Fortran 2008 mode, bounds/runtime checks, warnings, and strict
implicit-interface checking.

## Main modules

- `tsa_statistics`: moments, ACF/PACF/CCF, raw periodogram, harmonics, season
  indices, and exact runs calculations.
- `tsa_spectral`: R-compatible univariate and multivariate `spec.pgram` numerical
  paths (compact/full symmetric kernels, Daniell/modified-Daniell/Fejer/
  Dirichlet kernels, scalar or per-series taper, cross-spectra, coherence, and
  phase) plus Yule-Walker, Burg/Burg2, OLS, and ML `spec.ar` paths.
- `tsa_arma`: OLS AR fitting/selection, ARMA spectra, prewhitening, EACF,
  ARMA subset selection through the supplied `leaps` port, and Box-Cox AR
  likelihood profiling.
- `tsa_tests`: Ljung-Box/Box-Pierce, McLeod-Li, Keenan, Tsay, and AO/IO
  outlier diagnostics.
- `tsa_simulation`: QAR, GARCH, and SETAR/TAR simulation.
- `tsa_tar`: two-regime SETAR fitting by MAIC or CLS, optional MAIC regime
  order selection, skeleton/cycle calculations, simulation prediction, and
  Chan threshold likelihood-ratio inference.
- `tsa_arimax`: CSS/CSS-ML and diffuse state-space ML for ARIMA/seasonal
  ARIMA and ARIMAX, transfer/intervention regressors, R/TSA-style BFGS
  parameter scaling, direct SVD conditioning for all-free regression blocks,
  Gardner/AS154 stationary initialization, R-style `optimHess` covariance
  construction, simulation, and bootstrap support.
- `tsa_garch_diagnostics`: TSA's generalized Box diagnostic for fitted GARCH
  models using the supplied `tseries` port.
- `tsa`: convenience facade exporting the public API.

## Scope

Plotting, R formula/S3 dispatch, printed summaries, and other presentation-only
code are deliberately omitted. `docs/API_MAP.md` maps the R exports to the
Fortran API. `docs/TRANSLATION_STATUS.md` records the narrower exact-parity limitations
that remain rather than hiding them.

## Licensing

TSA is GPL (>= 2). The supplied tseries translation is GPL-2.0-only OR
GPL-3.0-only, so the combined bundle is distributable under GPL version 2 or
GPL version 3. See `LICENSE`, `LICENSE-GPL-2`, `LICENSE-GPL-3`, `NOTICE`, and
the license files retained inside `upstream/` and `vendor/`.
