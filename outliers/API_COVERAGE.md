# API coverage

Package status: substantial. 14 of 14 (100%) exported computational R functions are mapped.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `chisq.out.test` | `outliers_api` | `chisq_out_test` | substantial |
| `cochran.test` | `outliers_api` | `cochran_test`, `cochran_groups` | substantial |
| `dixon.test` | `outliers_api` | `dixon_test` | substantial |
| `grubbs.test` | `outliers_api` | `grubbs_test` | substantial |
| `outlier` | `outliers_api` | `outlier`, `outlier_mask` | substantial |
| `pcochran` | `outliers_distributions` | `pcochran` | complete |
| `pdixon` | `outliers_distributions` | `pdixon` | substantial |
| `pgrubbs` | `outliers_distributions` | `pgrubbs` | substantial |
| `qcochran` | `outliers_distributions` | `qcochran` | complete |
| `qdixon` | `outliers_distributions` | `qdixon` | substantial |
| `qgrubbs` | `outliers_distributions` | `qgrubbs` | substantial |
| `qtable` | `outliers_distributions` | `qtable` | complete |
| `rm.outlier` | `outliers_api` | `rm_outlier` | substantial |
| `scores` | `outliers_api` | `scores` | substantial |

## Mapping details

* `chisq.out.test`: NaN omission, selected endpoint, statistic and upper-tail probability; typed result replaces htest and labels.
* `cochran.test`: Variance/count, equal-size matrix, and integer-labelled grouped samples support outlying/inlying tests and ties; R formula evaluation and htest construction omitted.
* `dixon.test`: All six ratio types and automatic type, opposite and two-sided options; upstream tables and anomalies preserved; typed result replaces htest.
* `grubbs.test`: All three variants, G/U statistics, opposite and two-sided options; type 11 inverse uses algebra with the upstream clipping/bracket domain.
* `outlier`: Vector and matrix-column values plus logical masks; all selected ties retained; missing mask entries are false rather than R NA.
* `pcochran`: Scalar/elemental cumulative probability using the shared F distribution, with documented valid parameter domain.
* `pdixon`: Scalar/elemental reverse table interpolation for every Dixon ratio type, including original anomalous entries.
* `pgrubbs`: Scalar/elemental distribution values for all three variants; type 11 uses the exact algebraic inverse rather than R uniroot tolerance.
* `qcochran`: Scalar/elemental critical values using the shared F quantile, including probability endpoints.
* `qdixon`: Scalar/elemental forward and reverse table interpolation with original sample-size restrictions and source table anomalies.
* `qgrubbs`: Scalar/elemental critical values and reverse mode for types 10, 11, 20; table domains and type 11 clipping/bracket gap retained.
* `qtable`: Scalar/vector local cubic interpolation and endpoint linear extrapolation; paired sorting and tied-endpoint regression retained; local interior duplicate coordinates return NaN.
* `rm.outlier`: Vector and matrix-column removal or mean/median replacement, opposite choice and all ties; matrix removal returns typed ragged columns.
* `scores`: All five score types, probability, zero-probability and IQR-limit modes, vector/matrix columns; threshold results are real 0/1 with NaN for missing, no data-frame dispatch.

All R sources and manuals are preserved under upstream/. No exported computational functions were omitted from the denominator. Although qtable is called an internal routine in its manual, exportPattern exports it, so it is counted. R presentation and object construction are not counted separately.
