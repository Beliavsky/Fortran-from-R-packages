# R stats numerical parity reference

v0.5.0 uses the current R source tree as the reference for behavior TSA
delegates to `stats`, while the bundled TSA 1.3.1 source remains the reference
for TSA-specific transfer/intervention and diagnostic behavior.

Primary R source files consulted through v0.5.0:

- `src/library/stats/R/arima.R` -- xreg SVD rotation, regression starts/scales,
  ARIMA fit orchestration, and coefficient/covariance back-transforms.
- `src/library/stats/R/optim.R`, `src/library/stats/src/optim.c`, and
  `src/appl/optim.c` -- BFGS and `optimHess` numerical behavior.
- `src/library/stats/R/spectrum.R` -- `spec.ar`, `spec.pgram`, tapering,
  cross-periodograms, coherence, phase, df, and bandwidth calculations.
- `src/library/stats/R/ar.R` and `src/library/stats/src/burg.c` -- Yule-Walker,
  Burg/Burg2, OLS, and MLE AR order-selection conventions and Burg recursion
  used by `spec.ar`.
- `src/library/stats/R/kernel.R` -- compact symmetric `tskernel` coefficients,
  normalization, named kernels, df/bandwidth, and circular smoothing semantics.

The R source is GPL-2-or-later. No additional R source file is copied into the
Fortran package; translated behavior is implemented independently and covered
by regression tests.
