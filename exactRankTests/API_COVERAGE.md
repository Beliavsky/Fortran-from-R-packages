# Computational API coverage

## Translation coverage

Package status: **substantial**. **15 of 18 (83.33%)** exported computational R
functions have meaningful Fortran mappings. The fraction measures mapped
computational functions, not complete R compatibility. Generic entry points
and explicitly exported numerical S3 methods count separately, as required by
the translation convention. The three exported formula methods remain in the
denominator conservatively and are untranslated: `wilcox.exact.formula`,
`ansari.exact.formula`, `perm.test.formula`.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `wilcox.exact` | `exact_rank_tests` | `wilcox_exact` | partial |
| `wilcox.exact.default` | `exact_rank_tests` | `wilcox_exact` | partial |
| `ansari.exact` | `exact_rank_tests` | `ansari_exact` | partial |
| `ansari.exact.default` | `exact_rank_tests` | `ansari_exact` | partial |
| `perm.test` | `exact_rank_tests` | `perm_test` | partial |
| `perm.test.default` | `exact_rank_tests` | `perm_test` | partial |
| `dperm` | `exact_rank_distribution` | `permutation_distribution`, `simulate_distribution`, `dperm` | substantial |
| `pperm` | `exact_rank_distribution` | `permutation_distribution`, `simulate_distribution`, `pperm`, `dperm` | substantial |
| `qperm` | `exact_rank_distribution` | `permutation_distribution`, `simulate_distribution`, `qperm` | substantial |
| `rperm` | `exact_rank_distribution` | `rperm` | substantial |
| `cscores` | `exact_rank_scores` | `cscores`, `cscores_factor`, `cscores_surv` | substantial |
| `cscores.default` | `exact_rank_scores` | `cscores` | substantial |
| `cscores.factor` | `exact_rank_scores` | `cscores_factor` | substantial |
| `cscores.Surv` | `exact_rank_scores` | `cscores_surv` | substantial |
| `irank` | `exact_rank_scores` | `irank` | substantial |

The table and counts match `fpm.toml`. Detailed boundaries appear in
[API_COVERAGE.md](API_COVERAGE.md).

## Per-entry boundaries

- `wilcox.exact`: One/two-sample and paired exact/asymptotic tests, exact shift intervals and estimates. Asymptotic confidence intervals, formula processing and S3 objects are omitted. Undefined upstream estimates are normalized to NaN.
- `wilcox.exact.default`: One/two-sample and paired exact/asymptotic tests, exact shift intervals and estimates. Asymptotic confidence intervals, formula processing and S3 objects are omitted. Undefined upstream estimates are normalized to NaN.
- `ansari.exact`: Exact/asymptotic tests and exact scale intervals for finite nonzero data. Asymptotic confidence intervals, formula processing and S3 objects are omitted.
- `ansari.exact.default`: Exact/asymptotic tests and exact scale intervals for finite nonzero data. Asymptotic confidence intervals, formula processing and S3 objects are omitted.
- `perm.test`: One/two-sample and paired exact/asymptotic tests, including upstream noninteger score remapping. Permutation confidence sets, formula processing and S3 objects are omitted.
- `perm.test.default`: One/two-sample and paired exact/asymptotic tests, including upstream noninteger score remapping. Permutation confidence sets, formula processing and S3 objects are omitted.
- `dperm`: Density queried from a reusable exact or simulated distribution. Full support/mass arrays replace the R density data-frame option. Supplied uniforms replace R RNG state; resource limits return status.
- `pperm`: Inclusive one/two-sided tail probabilities and point mass via dperm; supplied uniforms replace R RNG state; resource limits return status.
- `qperm`: Left-continuous inverse CDF on a reusable exact/simulated distribution, with robust endpoint handling. Supplied uniforms replace R RNG state.
- `rperm`: Uniform fixed-size sampling without replacement through caller-supplied independent uniforms and Fisher-Yates selection; R RNG sequences are not reproduced.
- `cscores`: All seven numerical score families, factor codes, log-rank scores and integer scaling. Type-specific entry points replace S3 dispatch; nonfinite values are rejected.
- `cscores.default`: All seven numerical score families, factor codes, log-rank scores and integer scaling. Type-specific entry points replace S3 dispatch; nonfinite values are rejected.
- `cscores.factor`: Numeric scores for one-based factor codes 1 and 2; label/level construction is caller-managed.
- `cscores.Surv`: Log-rank scores including ties/censoring and optional integer scaling; explicit time/event arrays replace Surv objects. Undefined rescaling returns status.
- `irank`: Maximum integer ranks with validated optional one-based sorting order; finite data only.
