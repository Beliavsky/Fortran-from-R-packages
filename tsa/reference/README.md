# Cross-language parity fixtures

`generate_r_parity_v040.R` generates deterministic reference output from R's
`stats::arima` and, when installed, `TSA::arimax` for the regression-conditioning
cases added in v0.4.0. `generate_r_parity_v050.R` generates fixed-order and AIC
AR-method spectra plus custom-kernel/per-series-taper periodogram references
for the v0.5.0 spectral work.

The validation container used for this translation does not contain an R
runtime, so the script is included as a reproducible external parity harness
rather than embedding unverified R output in the repository.

The corresponding native Fortran assertions are in:

- `test/test_xreg_parity.f90`
- `test/test_optimizer_parity.f90`
- `test/test_linalg_parity.f90`
- `test/test_spectral_methods.f90`
