# Computational coverage

| Original R routine | Fortran routine | Status |
|---|---|---|
| `MOE` | `moe`, `fit_all_densities` | Numerical workflow translated; plotting and CSV side effects excluded |
| `approximate.max` | `approximate_max` | Translated |
| `bsm.objective` | `bsm_objective` | Translated |
| `compute.implied.volatility` | `compute_implied_volatility` | Translated with positive-volatility bracketing |
| `dew` | `dew` | Translated |
| `dgb` | `dgb` | Translated |
| `dmln` | `dmln` | Translated |
| `dmln.am` | `dmln_am` | Translated |
| `dshimko` | `dshimko` | Translated |
| `ew.objective` | `ew_objective` | Translated |
| `extract.am.density` | `extract_am_density` | Translated |
| `extract.bsm.density` | `extract_bsm_density` | Translated |
| `extract.ew.density` | `extract_ew_density` | Translated |
| `extract.gb.density` | `extract_gb_density` | Translated |
| `extract.mln.density` | `extract_mln_density` | Translated |
| `extract.rates` | `extract_rates` | Translated |
| `extract.shimko.density` | `extract_shimko_density` | Translated |
| `fit.implied.volatility.curve` | `fit_implied_volatility_curve` | Translated |
| `gb.objective` | `gb_objective` | Translated |
| `get.point.estimate` | `get_point_estimate` | Translated |
| `mln.am.objective` | `mln_am_objective` | Translated with corrected indexing |
| `mln.objective` | `mln_objective` | Translated |
| `pgb` | `pgb` | Translated |
| `price.am.option` | `price_am_option` | Translated |
| `price.bsm.option` | `price_bsm_option` | Translated |
| `price.ew.option` | `price_ew_option` | Translated |
| `price.gb.option` | `price_gb_option` | Translated |
| `price.mln.option` | `price_mln_option` | Translated |
| `price.shimko.option` | `price_shimko_option` | Translated using the equivalent closed form |

## Excluded infrastructure

- R plotting, graphics devices, and PDF reports
- CSV-writing side effects from `MOE`
- `lm`, `summary.lm`, formula, list, and data-frame presentation objects
- R package namespace machinery
- Compiled use of the four bundled `.rda` data sets

The supplied R package, including data and documentation, is retained in
`original/RND-1.2` for provenance and comparison.
