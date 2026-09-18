# API coverage

Coverage basis: distinct exported computational R functions and exported computational S3 methods in the upstream
`NAMESPACE`. Plotting and printing methods are excluded. Internal helpers such as `.ECM`, `.GssERC`, `.EqNPC`,
`.EqERC`, weighted-ECDF helpers, and row/assignment utilities are not counted separately.

Package status: **substantial**.

**12 of 12 (100.0%)** exported computational functions/methods have meaningful Fortran mappings. This percentage
measures mapped computational functions, not complete R interface compatibility or bit-for-bit reproduction.

| R function | Fortran module | Fortran procedure(s) | Status |
| --- | --- | --- | --- |
| `InitClust` | `otrimle_mod` | `init_clust` | substantial |
| `otrimle` | `otrimle_mod` | `otrimle_fit_grid` | substantial |
| `rimle` | `otrimle_mod` | `rimle` | substantial |
| `kerndensmeasure` | `otrimle_mod` | `kerndensmeasure` | substantial |
| `kmeanfun` | `otrimle_mod` | `kmeanfun` | complete |
| `ksdfun` | `otrimle_mod` | `ksdfun` | complete |
| `kerndensp` | `otrimle_mod` | `kerndensp` | substantial |
| `kerndenscluster` | `otrimle_mod` | `kerndenscluster` | substantial |
| `generator.otrimle` | `otrimle_mod` | `generator_otrimle` | substantial |
| `otrimleg` | `otrimle_mod` | `otrimleg` | substantial |
| `otrimlesimg` | `otrimle_mod` | `otrimlesimg` | substantial |
| `summary.otrimlesimgdens` | `otrimle_mod` | `summarize_otrimlesimgdens` | substantial |

Material differences are documented in `README.md`, `NOTICE.md`, and the corresponding `fpm.toml` mapping notes.
