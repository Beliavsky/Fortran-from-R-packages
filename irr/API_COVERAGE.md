# Computational API coverage

Status: substantial. 19 of 19 (100%) distinct exported computational functions.
The inventory follows `exportPattern(".")` in the original NAMESPACE and all
21 top-level R functions; the two print methods are outside the denominator.
Mappings match fpm.toml and README.md. All listed procedures are public through
both their implementation module and the `irr` facade.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `agree` | `irr_agreement` | `agree` | substantial |
| `bhapkar` | `irr_marginal` | `bhapkar`, `bhapkar_table` | substantial |
| `finn` | `irr_agreement` | `finn` | substantial |
| `icc` | `irr_icc` | `icc` | substantial |
| `iota` | `irr_agreement` | `iota` | substantial |
| `kappa2` | `irr_kappa` | `kappa2`, `kappa2_table` | substantial |
| `kappam.fleiss` | `irr_kappa` | `kappam_fleiss` | substantial |
| `kappam.light` | `irr_kappa` | `kappam_light` | substantial |
| `kendall` | `irr_agreement` | `kendall` | substantial |
| `kripp.alpha` | `irr_agreement` | `kripp_alpha` | substantial |
| `maxwell` | `irr_agreement` | `maxwell` | substantial |
| `meancor` | `irr_agreement` | `meancor` | substantial |
| `meanrho` | `irr_agreement` | `meanrho` | substantial |
| `N.cohen.kappa` | `irr_designs` | `n_cohen_kappa` | substantial |
| `N2.cohen.kappa` | `irr_designs` | `n2_cohen_kappa` | substantial |
| `rater.bias` | `irr_marginal` | `rater_bias`, `rater_bias_ratings` | substantial |
| `relInterIntra` | `irr_designs` | `rel_inter_intra` | substantial |
| `robinson` | `irr_agreement` | `robinson` | substantial |
| `stuart.maxwell.mh` | `irr_marginal` | `stuart_maxwell_mh`, `stuart_maxwell_ratings` | substantial |

## Mapping notes

- **agree**: Percentage agreement with numeric tolerance and complete-case deletion.
- **bhapkar**: Bhapkar marginal-homogeneity statistic and chi-square tail; explicit score and count-table entry points.
- **finn**: Both ANOVA models, coefficient, infinite-numerator-df F statistic and p-value.
- **icc**: One-/two-way, consistency/agreement, single/average, null reliability, F inference and confidence bounds; Fortran keyword kind replaces R type.
- **iota**: Quantitative or nominal multivariate arrays and optional standardization; explicit cross-variable listwise NaN deletion.
- **kappa2**: Unweighted, linear, squared and arbitrary R-normalized numeric weights; numeric categories, no factor metadata.
- **kappam.fleiss**: Fleiss and exact Conger estimates with unrounded category details; R detail-table rounding omitted.
- **kappam.light**: Mean pairwise kappas and upstream variance formula; negative variance retains coefficient and returns degenerate inference.
- **kendall**: Average ranks, optional tie correction, chi-square statistic and tail.
- **kripp.alpha**: All four distances, pairwise missingness and coincidence matrix. R complete-data weights default; optional corrected normalization.
- **maxwell**: Binary coefficient with explicit category count validation.
- **meancor**: Mean Pearson correlations, Fisher transform and upstream inference; dropped perfect-pair count.
- **meanrho**: Mean Spearman correlations with average ranks, Fisher transform and upstream inference.
- **N.cohen.kappa**: Binary sample size, both test tails and optional raw size; infeasible probability tables rejected.
- **N2.cohen.kappa**: Continuous LP variance maximization for 2-10 categories via sibling lpSolve; solver failures return status.
- **rater.bias**: Ordered disagreement ratio and chi-square test; explicit table/rating entry points replace R shape guessing.
- **relInterIntra**: All numeric reliability, F, p-value, lower-bound and SEM vectors for balanced complete data; R formula/labels omitted.
- **robinson**: Robinson A based on two-way sums of squares.
- **stuart.maxwell.mh**: Upstream equal-marginal removal default; optional retention of discordant equal-marginal categories; singular covariance returns status.

## Repeated-measure vector order

For R raters and M repeats, `rel_inter_intra_result` contains:

- `rohat(4+2*R)`: inter random, intra random, inter fixed, intra fixed,
  then R per-rater intra random and R per-rater intra fixed values.
- `fs(2+R)` and `pvalue(2+R)`: inter, pooled intra, then per-rater intra.
- `lowvalue(3+R)`: inter random, inter fixed, pooled intra, then per-rater intra.
- `sem(3+R)`: pooled intra, R per-rater intra, inter random, inter fixed.
- Extra diagnostics: `mean_squares(4+R)` holds MSS, MSR, MSSR, MSE and per-rater
  MSE; `variance_random(4)` and `variance_fixed(4)` hold subject, rater,
  interaction, residual components; `nu(2)` holds inter-bound denominator df.

Negative estimated variance components are retained as upstream ANOVA estimates.
Undefined standard errors or confidence bounds are NaNs; consult status.

## Omitted interface-only operations

`print.irrlist`, `print.icclist`, textual diagnostics, labels, data-frame/factor
construction and coercion are not translated. Datasets are not computational
functions. No exported computational function is left untranslated.
