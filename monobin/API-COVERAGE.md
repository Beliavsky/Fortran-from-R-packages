# Computational API coverage

Coverage basis: distinct exported computational R functions in monobin 0.2.4
NAMESPACE. Total 7; translated 7; fraction 1.0; untranslated: none. Package and
all seven mappings are **substantial**, reflecting material R interface changes.

| R function | Source | Public Fortran procedure | Numerical operation |
|---|---|---|---|
| pct.bin | upstream/R/01_PCT_BINNING.R | monobin::pct_bin | Iterative percentile discretization and trend check |
| iso.bin | upstream/R/02_ISO_BINNING.R | monobin::iso_bin | Isotonic regression plus count/rate corrections |
| woe.bin | upstream/R/03_WOE_BINNING.R | monobin::woe_bin | Isotonic bins plus adjacent WoE-gap correction |
| cum.bin | upstream/R/04_CUM_BINNING.R | monobin::cum_bin | Cumulative response-mean pooling |
| sts.bin | upstream/R/05_STS_BINNING.R | monobin::sts_bin | Adjacent significance-test merging |
| ndr.bin | upstream/R/06_NDR_BINNING.R | monobin::ndr_bin | Backward nested-dummy regression cut selection |
| mdt.bin | upstream/R/07_MDT_BINNING.R | monobin::mdt_bin | Monotone decision tree |

All public procedures above are in `src/monobin.f90`. Public `fit_bins` and
`apply_bins` are typed dispatch/reuse extensions and add no R coverage count.
Private helpers replace dplyr group operations with sufficient-statistic arrays.
Internal R helpers, printed labels, plotting, data-frame construction, and the
`gcd` data object are not additional exported computational functions. No fitting
function is excluded because its upstream interface uses data frames.

README's translation section and the manifest carry matching counts and mapping
status. Source details, optional controls, missing/special data handling,
infinite evidence, regression edge cases, and preserved Hmisc/MDT endpoint
quirks are documented there rather than claimed as exact R object compatibility.
