# Validation

The v0.5.0 package was validated on 2026-08-18 with GNU Fortran 14.2.0. FPM was
not installed in the validation environment, so the FPM source tree was
compiled directly in dependency order from the standalone sources under
`src/`.

Compiler checks used:

```text
-std=f2008 -O0 -g -Wall -Wextra -Wimplicit-interface
-fcheck=all -ffpe-trap=invalid,zero,overflow
```

Clean-build results:

```text
test_arimax_features: PASS
test_arimax_parity: PASS
test_bootstrap_parity: PASS
test_diffuse_ml_parity: PASS
test_linalg_parity: PASS
test_optimizer_parity: PASS
test_spectral_methods: PASS
test_spectral_parity: PASS
test_tar_parity: PASS
test_tsa: PASS
test_xreg_parity: PASS
demo_tsa.f90: PASS
```

The v0.5.0 additions are specifically covered by:

- fixed-order Yule-Walker, Burg/Burg2, OLS, and ML AR spectra;
- method-specific AIC order-selection paths;
- compact `tskernel` expansion and named-kernel normalization;
- equality of compact and explicitly expanded smoothing weights;
- multivariate per-series taper autospectra, df corrections, and bounded
  coherence;
- the v0.4 SVD/xreg, Hessian, covariance, and linear-algebra parity cases.

The complete suite also continues to cover transfer/intervention ARIMAX,
diffuse integrated/seasonal ML, bootstrap parity, TAR, spectra, diagnostics,
and the original TSA numerical tests.

All Fortran files under `src/`, `test/`, and `example/` are free-format.
