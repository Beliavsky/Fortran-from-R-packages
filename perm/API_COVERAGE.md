# Computational API coverage

The denominator is all **20 distinct exported computational R functions** in
perm 1.0-0.4, as listed by the official NAMESPACE. All 20 are mapped. Package
status is **substantial**, not a claim of complete R compatibility.

The exported `permControl` is included because it performs input validation and
constructs computational controls. `getcnt` is included because it computes
variable-length network edge-index concatenations; it is not a presentation-only
accessor. The three `perm*` S3 generics count once each: their numerical default
workflow is translated. The registered default/formula methods are not separately
exported, so they do not increase this denominator. `print.mchtest` is unexported
and presentation-only. No datasets or internal helpers increase the count.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `chooseMatrix` | `perm_api` | `choose_matrix` | complete |
| `permTS` | `perm_api` | `perm_ts` | substantial |
| `permKS` | `perm_api` | `perm_ks` | substantial |
| `permTREND` | `perm_api` | `perm_trend` | substantial |
| `getcnt` | `perm_api` | `getcnt` | complete |
| `methodRuleTS1` | `perm_api` | `method_rule_ts1` | substantial |
| `methodRuleKS1` | `perm_api` | `method_rule_ks1` | substantial |
| `methodRuleTREND1` | `perm_api` | `method_rule_trend1` | substantial |
| `twosample.exact.ce` | `perm_api` | `twosample_exact_ce` | substantial |
| `twosample.exact.mc` | `perm_api` | `twosample_exact_mc` | substantial |
| `twosample.exact.network` | `perm_api` | `twosample_exact_network` | substantial |
| `twosample.pclt` | `perm_api` | `twosample_pclt` | substantial |
| `ksample.pclt` | `perm_api` | `ksample_pclt` | substantial |
| `ksample.exact.mc` | `perm_api` | `ksample_exact_mc` | substantial |
| `trend.pclt` | `perm_api` | `trend_pclt` | substantial |
| `trend.exact.mc` | `perm_api` | `trend_exact_mc` | substantial |
| `calcPvals` | `perm_api` | `calc_pvals` | substantial |
| `permControl` | `perm_api` | `perm_control` | substantial |
| `calcPvalsMC` | `perm_api` | `calc_pvals_mc` | substantial |
| `pCI` | `perm_api` | `p_ci` | substantial |

## Per-function details

- **chooseMatrix** (`upstream/R/chooseMatrix.R`): Enumerates zero/one selection rows in upstream order, with explicit allocation and combinatorial limits.
- **permTS** (`upstream/R/permTS.R`): Numerical two-sample workflow and method selection; formulas, S3 objects, R warnings and custom method-rule callbacks are omitted.
- **permKS** (`upstream/R/permKS.R`): Numerical k-sample workflow; integer category labels replace R factors. This is not a Kolmogorov-Smirnov test; formulas, S3 objects and callbacks are omitted.
- **permTREND** (`upstream/R/permTREND.R`): Numerical trend workflow and correlation estimate; formulas, S3 objects, R warnings and callbacks are omitted.
- **getcnt** (`upstream/R/getcnt.R`): Concatenates outgoing network-edge indices, with one-based indices and explicit shape and bounds validation.
- **methodRuleTS1** (`upstream/R/methodRules.R`): Distinct-score-count method thresholds; exact=-1/0/1 replaces NULL/FALSE/TRUE. One-distinct-score groups use the first bound rather than the upstream index-zero error.
- **methodRuleKS1** (`upstream/R/methodRules.R`): Smallest-group method threshold; integer labels replace R factors and exact=-1/0/1 replaces NULL/FALSE/TRUE.
- **methodRuleTREND1** (`upstream/R/methodRules.R`): Sample-size method threshold, with exact=-1/0/1 replacing NULL/FALSE/TRUE and explicit finite-data validation.
- **twosample.exact.ce** (`upstream/R/algorithms.R`): Complete binary-label enumeration or supplied selection rows, with both two-sided conventions, significant-digit comparisons and resource guards.
- **twosample.exact.mc** (`upstream/R/algorithms.R`): Fixed-count permutations with add-one tails and Monte Carlo intervals; local Park-Miller RNG differs from R and requires an explicit stream for continuation.
- **twosample.exact.network** (`upstream/R/algorithms.R`): Exact grouped-score occupancy paths weighted by binomial multiplicities, using log weights to avoid overflow. Arbitrary real scores and ties are retained; resource limits are explicit.
- **twosample.pclt** (`upstream/R/algorithms.R`): Permutation normal approximation, with the larger integer label selected; stable upper tails and explicit degenerate-data status replace cancellation or NaN.
- **ksample.pclt** (`upstream/R/algorithms.R`): Centered between-group quadratic statistic and chi-square approximation; integer labels replace factors and constant scores return a degenerate status.
- **ksample.exact.mc** (`upstream/R/algorithms.R`): Fixed-count k-sample permutations and upper-tail uncertainty; local RNG and power-of-two scaling differ from R.
- **trend.pclt** (`upstream/R/algorithms.R`): Correlation-based permutation normal approximation with stable scaling and explicit degenerate-data status.
- **trend.exact.mc** (`upstream/R/algorithms.R`): Fixed-count linear-statistic permutations and both two-sided conventions; local RNG, scaled statistics and explicit continuation replace R global RNG handling.
- **calcPvals** (`upstream/R/algorithms.R`): Weighted centered-statistic tails with ties-to-even significant-digit comparisons. The unused R denom argument is omitted; extreme finite scales are handled safely.
- **permControl** (`upstream/R/permControl.R`): Constructs and validates sampling, rounding, confidence, supplied-matrix and tail-convention controls, with additional explicit resource bounds.
- **calcPvalsMC** (`upstream/R/algorithms.R`): Add-one Monte Carlo tails and Clopper-Pearson intervals; both doubled central interval endpoints are capped at one to repair reversed all-ties intervals upstream.
- **pCI** (`upstream/R/algorithms.R`): Equal-tailed Clopper-Pearson interval through shared beta/F numerics. Stable tail symmetry avoids upstream probability-rounding loss at confidence near one; trial count is bounded.

## Interface differences that apply throughout

No formula/data-frame evaluation, missing-data recycling, S3 dispatch, printed
objects, R warnings, character/factor conversion, implicit global RNG state or
custom method-rule callbacks are implemented. Supply finite arrays with matching
shapes and integer category labels. Select a custom method outside the library
and pass its explicit name. `permKS` means **k-sample**.

Status codes and all practical limits are listed in README.md. Limit statuses are
failures, never approximate p-values disguised as exact results. Monte Carlo
procedures run a fixed number of samples and have no sequential stopping rule.
