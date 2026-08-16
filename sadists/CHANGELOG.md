# Changelog

## 0.2.0 - 2026-08-15

- Refactor the approximation layer to depend on standalone
  `pdqutils-fortran` 0.1.0 instead of duplicating PDQutils algorithms.
- Replace the former Edgeworth/Cornish-Fisher implementation in
  `sadists_approximations` with a thin compatibility adapter.
- Reuse PDQutils normal PDF/CDF/quantile support from `sadists_special`.
- Reuse PDQutils moment-to-cumulant and cumulant-to-moment implementations
  through API-compatible sadists wrappers.
- Vendor `pdqutils-fortran` and add it as a local FPM dependency so the source
  archive remains self-contained.
- Add `test_pdqutils_integration` to verify density, CDF, quantile, and
  moment/cumulant agreement across the dependency boundary.
- Preserve the v0.1.0 public sadists API and numerical results.

## 0.1.0 - 2026-08-15

- Initial modern Fortran/FPM translation of sadists 0.2.6 computational code.
- Port all 12 exported distribution families and 48 d/p/q/r operations.
- Translate required Edgeworth and Cornish-Fisher/AS269 machinery.
- Replace hypergeo and orthopolynom runtime dependencies with standalone
  moment algorithms.
- Add adaptive noncentral log-chi-square Poisson-mixture moments.
- Add standalone normal, gamma, Poisson, and noncentral chi-square RNGs.
- Add stable upper/log-tail normal quantiles and support endpoint handling.
- Preserve R-style recycling among multi-component parameter vectors.
- Omit the Shiny-only `runExample` function.
- Add regression, transformation, numerical, and RNG tests plus an example.
