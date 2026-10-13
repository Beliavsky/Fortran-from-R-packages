# API coverage

Status: **substantial**. Coverage basis: **exported computational R
functions**. **2 of 2 (100%)** functions are mapped, as recorded in
`fpm.toml`. This count does not assert complete R compatibility.

| R function | Source | Public Fortran procedure | Status |
| --- | --- | --- | --- |
| `pr.curve` | `upstream/R/PRROC.R` | `prroc_api:pr_curve` in `src/prroc_api.f90` | substantial |
| `roc.curve` | `upstream/R/PRROC.R` | `prroc_api:roc_curve` in `src/prroc_api.f90` | substantial |

Both mappings include weighted and unweighted observations, binary/soft
labels, sorted input, tied thresholds, optional curve coordinates, and
maximum/minimum/random references. PR includes exact nonlinear integrated
area, optional hard-label Davis–Goadrich area, and interpolated curve points.
ROC includes weighted trapezoidal area.

No exported computational function is untranslated. `ggprroc` is an exported
plotting/presentation helper, excluded from the coverage basis. Registered
`plot.PRROC` and `print.PRROC` methods, internal helpers (`compute.pr`,
`compute.roc`, `check`, and curve-allocation helpers), and S3 list/class
handling are not independently counted.

The Fortran result uses a common `auc` component rather than R's separate
`auc.integral`/`auc` names. See README for all arguments and types. Private
helpers are implementation details and are not claimed as mapped public APIs.

Material differences: explicit validation and statuses; no R recycling or
presentation interface; removal of zero-mass thresholds and random-matrix
NA padding; mathematically correct random PR prevalence for unequal-length
weights; finite-weight requirement; bounded curve allocation; stable
class-wise normalization and nonlinear integration. R lazy score/weight
default evaluation order is retained. The reference extrema retain upstream's
ordering constructions, including their limitations for independent,
non-complementary weights on matching score vectors.

Validation comprises actual-R generated cases and independent mathematical
properties. Upstream curve rows changed by the documented zero-weight/padding
policy are not falsely claimed as exact row-for-row matches.
