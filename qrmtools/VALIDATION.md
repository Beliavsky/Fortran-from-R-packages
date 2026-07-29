# Validation report

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Checked build: `-O0`, bounds and runtime checks, backtraces, strict interface
  and conversion diagnostics, and warnings treated as errors
- Optimized build: `-O2` with the same diagnostics and runtime checks

FPM was not installed in the validation container. The manifest was parsed as
TOML and the exact source graph was compiled directly in FPM dependency order.

## Test suites

```text
test_bounds_allocation: PASS
test_brownian_black_scholes: PASS
test_diagnostics_garch: PASS
test_distributions_evt: PASS
test_geometric_composite: PASS
```

The demo and both examples also compile and run under checked and optimized
builds.

## Numerical references

The suites include fixed references for:

- GEV and GPD densities, CDFs, quantiles, and expected shortfall;
- Pareto and GPD-tail probability identities;
- GPD moment fitting and GEV/GPD maximum-likelihood smoke tests;
- Hill estimates and tied-observation mean excess;
- standardized Student-t VaR and ES;
- Black-Scholes calls, puts, and Greeks;
- deterministic Brownian increments and exact de-Browning recovery;
- return transformations and inversion;
- hierarchical matrices and risk allocation;
- column-preserving rearrangement and adaptive bound construction;
- Mahalanobis distance identities and multivariate diagnostic ranges;
- a fixed reparameterized GARCH log-likelihood;
- composite-distribution CDF/quantile inversion;
- geometric VaR and expectile objective improvement.

Reference values were calculated independently with NumPy/SciPy or direct
closed-form calculations.

## Release audit

The release validation checks:

- valid `fpm.toml` syntax;
- ASCII-only translated source and documentation;
- free-form source lines no longer than 132 columns;
- `implicit none` in every translated Fortran program unit;
- `GPL-3.0-or-later` SPDX headers in every translated Fortran file;
- absence of build products in the release tree;
- original, translated, and archive SHA-256 manifests;
- a clean rebuild and rerun from the exact final ZIP archive.
