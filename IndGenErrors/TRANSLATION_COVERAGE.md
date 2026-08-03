# Translation coverage

## Exported computational routines

| R export | Fortran routine | Status |
|---|---|---|
| `cvm_2series` | `cvm_2series` | translated |
| `cvm_3series` | `cvm_3series` | translated |
| `crosscor_2series` | `crosscor_2series` | translated |
| `crosscor_3series` | `crosscor_3series` | translated |
| `crossdep_2series` | `crossdep_2series` | translated |
| `crossdep_3series` | `crossdep_3series` | translated |

## Plotting-only exports

| R export | Status |
|---|---|
| `dependogram` | omitted |
| `CrossCorrelogram` | omitted |

## Internal computational coverage

- maximum-rank construction with ties
- circular two- and three-series lag kernels
- pair and triple Cramer-von Mises Mobius statistics
- finite-sample CDF tables for `n=50` and `n=100`
- bias-corrected W statistics
- Fisher p-value combinations
- cumulant and Edgeworth approximations
- chi-square tail probabilities
- Spearman, van der Waerden, and Savage Mobius scores
