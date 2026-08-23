# Translation coverage

## Implemented directly in Rfast2-fortran

### Array and utility kernels
`Quantile`, `rowQuantile`, `colQuantile`, `trim.mean`, `rowTrimMean`,
`colTrimMean`, `Intersect`, `Merge`, `colGroup`, `jack.mean`,
`coljack.means`, `rowjack.means`, `colmeansvars`, triangular checks,
`lud`, column accuracy/sensitivity/specificity/precision/F-score/FMI/FB-score,
MSE/MAE and KL-style column metrics.

### Random and sampling
`Runif`, `Sample.int`, `Sample`, `rbeta1`, plus the underlying distribution
RNG kernels used by Rfast2.

### Tests/statistics
`jbtest`, `jbtests`, `cor_test`, `covar`, `pooled.colVars`,
`empirical.entropy`, `pinar1`, `colpinar1`, `circ.cor1`, `circ.cors1`, `km`,
`moranI`, `wald.poisrat`, `walter.ci`, `wlsmeta`, `refmeta`, `silhouette`,
`perm.ttest1`, `perm.ttest2`, `boot.ttest1`, and the univariate core of
`eqdist.etest`/`normal.etest`.

### MLE
`gammapois.mle`, `halfcauchy.mle`, `cauchy0.mle`, `kumar.mle`,
`powerlaw.mle`, `zigamma.mle`, `zil.mle`, `ziweibull.mle`, `simplex.mle`,
`gnormal0.mle`, `unitweibull.mle`, `cbern.mle`, `sp.mle`; column-wise beta,
Cauchy, half-Cauchy, half-normal, lognormal, logit-normal, Borel, power-law,
unit-Weibull and SP fits.

### Regression/model fitting
`cls`, `covrob.lm`, the type-1 numerical path of `het.lmfit`, `gammareg`,
`gammaregs`, `weib.regs`, `multinom.reg`, grouped `binom.reg`, `ztp.reg`,
`sp.logiregs`, `batch.logistic`, `tobit.reg`, and reusable logistic/Poisson
wrappers.

### Multivariate
`pca`, `pcr`, `depth.mahala`, `leverage`, `diffic`, `discrim`, and a
covariance-distance kernel.

## Supplied by the active Rfast-fortran dependency

Rfast2 imports many operations directly from Rfast. The vendored dependency
provides modern Fortran implementations for the corresponding core domains:
array reductions/order statistics, matrix algebra, distances and energy
statistics, common distributions and MLEs, GLMs and regression, t/F/chi-square
and categorical tests, repeated-measures algorithms, directional statistics,
naive Bayes, score tests, empirical likelihood, OMP/selection primitives,
PC skeleton search, k-nearest-neighbor primitives, and multivariate MLEs.

These retain their Rfast-fortran names rather than being duplicated solely to
mirror an Rfast2 R wrapper.

## Remaining worthwhile numerical targets

The main Rfast2-specific numerical algorithms not yet ported one-for-one are:

- full FBED, MMPC/MMPC2, FEDHC/MMHC wrapper semantics and hash caches
- `pc.sel` and the complete multi-family `omp2` control flow
- censored Poisson and censored-Weibull MLE/regression
- Purkayastha and truncated-Cauchy MLEs
- the SPML-specific MLE/regression family and related naive Bayes wrappers
- GEE and clustered-regression orchestration
- full circular ANOVA family (`hcf`, `het`, `lr`, `embed`) and multi-sample
  directional fitting wrappers
- full James/Hotelling/resampling wrapper family
- all specialized naive-Bayes families and CV wrappers that are not already
  represented by Rfast-fortran's generic NB modules
- `frechet.nn`, `big.knn`/`bigknn.cv` and Rnanoflann-backed wrapper semantics
- specialized covariance equality/likelihood tests and RIAG
- full `lm.boot`, `wild.boot`, `lm.drop1`, `lm.bsreg` wrapper behavior
- `reg.mle.lda`, `fisher.da` and their CV wrappers

Plotting, benchmarking output, R environments/hash objects, formula parsing,
data frames, and R-level parallel orchestration are not numerical targets.
