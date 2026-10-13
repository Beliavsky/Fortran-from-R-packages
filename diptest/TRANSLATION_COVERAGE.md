# Computational export inventory

Source: diptest 0.77-2 `upstream/NAMESPACE` and its R definitions.
Basis: distinct exported computational R functions.
Package status: **substantial**. **2 of 2 (100%)** are translated substantially;
none are untranslated. This is an operation-coverage count, not full R-interface
or bitwise numerical compatibility.

| R export | Source | Public Fortran implementations | Status |
|---|---|---|---|
| `dip` | `upstream/R/dip.R` | `diptest_api: dip`, `dip_compute` | substantial |
| `dip.test` | `upstream/R/dipTest.R` | `diptest_api: dip_test`, `dip_test_simulate`, `dip_p_value` | substantial |

`dip` covers the statistic and both levels of GCM/LCM diagnostics, including
`min.is.0`. Debug output, R calls/classes and string-valued `full.result` dispatch
are replaced by a typed result with all diagnostics. Degenerate diagnostics and
extreme-range arithmetic are explicitly documented extensions.

`dip.test` covers NaN omission, the original simulation quantile table,
interpolation, small-n handling, and simulated exceedance proportions. R's RNG
is not translated: supply exact null samples for reproducible comparisons or use
the intrinsic Fortran RNG. `htest` construction, printed warnings/messages and
expression-derived data names are omitted.

## Exclusions and helper accounting

- `print.dip` and `plot.dip` are registered, unexported presentation methods and
  are excluded from the denominator.
- `getCM` is unexported, and its full-hull reconstruction is included in
  `dip_compute`; it is not counted as a separate R function.
- `aLine`, `rdRDS`, datasets and debug/registration interfaces are not exported
  computational functions. Table loading is replaced by compile-time data.
- Imported `graphics` functions are plotting-only. The needed sorting,
  NaN-removal, interpolation and RNG behavior of the imported `stats` helpers is
  implemented or reused as documented; the stats package is not vendored.

The two mapping entries and integer/fraction totals agree with `fpm.toml` and
README.md.
