# Computational API inventory

Coverage basis: **exported computational R functions** in HDInterval 0.2.4 `NAMESPACE`.
Package status: **substantial**, 2 of 2 mapped (100%), no unmapped exports.
This counts function mappings and does not assert full R compatibility.

| R export | Source | Public Fortran implementation | Status |
|---|---|---|---|
| `hdi` | `upstream/R/hdi.R`; methods in `hdiVector.R`, `hdi_density.R`, `hdi_function.R` | `hdinterval`: `hdi_samples`, `hdi_columns`, `hdi_density`, `hdi_function` | substantial |
| `inverseCDF` | `upstream/R/inverseCDF.R` | `hdinterval`: generic `inverse_cdf`, `inverse_cdf_scalar`, `inverse_cdf_vector` | substantial |

## Registered S3 methods (not additional exported functions)

| Registered method | Numerical coverage or adaptation |
|---|---|
| `hdi.default` | Numeric draws via `hdi_samples`; nonnumeric R objects are outside the typed API. |
| `hdi.matrix` | Independent columns via `hdi_columns`. |
| `hdi.data.frame` | Numeric columns via `hdi_columns`; heterogeneous R columns excluded. |
| `hdi.list` | Call an appropriate Fortran routine for each element; heterogeneous recursive dispatch excluded. |
| `hdi.density` | Threshold, ties, split modal intervals and contiguous fallback via `hdi_density`. |
| `hdi.function` | Quantile callback via `hdi_function`; local golden-section search replaces Brent. |
| `hdi.mcarray` | Caller reshapes/permutates the array to draw-by-variable columns. |
| `hdi.mcmc`, `hdi.mcmc.list` | Caller extracts and combines draw matrices; no `coda` object machinery. |
| `hdi.bugs`, `hdi.rjags`, `hdi.runjags`, `hdi.jagsUI` | Caller extracts samples; package object handling excluded. |

`hdiVector` and `checkCredMass` are unexported helpers. Their numerical operations
are included without counting them as exports. All six upstream R source files
are retained unmodified for audit and actual-R reference regeneration.

## Deliberate numerical and interface differences

- One named real kind, pure typed procedures and explicit status codes replace
  R objects, S3 dispatch, attributes, result names, exceptions and warnings.
- Samples preserve the exact `floor(n * mass)` spacing, first-width tie policy,
  NaN removal and infinity fallback of the R helper.
- Density grids require finite nonnegative ordinates and regular coordinates.
  Uniform-grid ordinate sums, threshold ties and strict cumulative comparisons
  match the upstream method. Scaling prevents overflow but may change a cutoff
  decision exactly at floating-point rounding boundaries.
- Nonsplit density fallback sets `is_hdi=.false.` rather than emitting a warning.
  A coarse grid with no candidate window returns a failure status and an empty
  result rather than NaN endpoints accompanied by R warnings.
- Callback dummy arguments are associated procedure pointers with `INTENT(IN)`.
  Parameters can be captured by pure host-associated wrappers instead of `...`.
- Quantile minimization is golden-section rather than R Brent; both depend on
  a unimodal width objective for a global optimum. Flat-width tie locations are
  not promised to match. Default tail tolerance remains `1e-8`.
- CDF inversion preserves `[1e-6,1-1e-6]` clipping, default `[-10,10]` starting
  bracket and the default `epsilon**0.25` tolerance. Symmetric bracket growth and
  bisection replace R's extension/Brent path. Optional budgets, bounds and
  per-element status codes are added. Bad vector elements do not discard good
  ones. `log.p` and `lower.tail` remain unavailable, as in upstream.

The README documents input domains, result shapes, every status code and
accuracy checks. The manifest and README coverage counts agree with this file.
