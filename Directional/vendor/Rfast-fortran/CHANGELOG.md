# Changelog

## 0.3.0

- Preserved the complete v0.1.0 and v0.2.0 public Fortran APIs.
- Completed all upstream `omp` response-family branches, including quasi families, normal-log, Gamma, Weibull, multivariate and multinomial selection.
- Corrected v0.2 OMP stopping semantics so the first selected model is tested against the upstream threshold before another variable is added.
- Added random-intercept MLE/regression, optional BLUP-style random effects, balanced `rint.regbx`, column wrappers, and batched `rint.regs`.
- Added normal-log and Weibull regressions, Tobit MLE, and ordinal logit/probit MLE.
- Added batched Gamma, inverse-Gaussian, quasi-Poisson, proportion, multinomial, normal-log and Poisson-family regression testing helpers.
- Added a native Pearson/Spearman PC-skeleton implementation with conditioning/separation sets.
- Added four v0.3 test suites and `example/v03_features.f90`.
- Expanded documented NAMESPACE coverage from 240 to 257 ported/consolidated exports; `pc.skel` is the sole partial export.

## 0.2.0

- Preserved the complete v0.1.0 public Fortran API.
- Added all ten Rfast score-test families and batched score calculations.
- Added multivariate normal/lognormal/t, Dirichlet, and inverse-Dirichlet MLEs.
- Added one-/two-sample empirical likelihood and exponential empirical likelihood, including multivariate EEL.
- Added repeated-measures ANOVA/line tests and variance-component MOM estimators.
- Added balanced and unbalanced variance-component MLEs with optional random effects; the unbalanced path translates upstream `gold_rat3`.
- Added Gamma, inverse-Gaussian, quasi-Poisson, proportion, multinomial, and spatial-median regressions.
- Added OMP/BIC linear selection, partial-correlation BIC forward selection, and logistic/Poisson GLM OMP.
- Added four v0.2 regression suites covering score tests, multivariate/EL methods, repeated measures/variance components, and regression/selection.
- Added `example/v02_features.f90`.
- Expanded documented NAMESPACE coverage from 204 to 240 ported/consolidated exports, with one additional partial export (`omp`).

## 0.1.0

- Initial modern Fortran 2018/FPM computational-core port of Rfast 2.1.5.2.
- Added row/column descriptive and ordering kernels, special functions, matrix algebra and distance/energy statistics.
- Added multivariate density/RNG support using the vendored zigg-fortran RNG.
- Added broad scalar and column-wise MLE coverage.
- Added linear/logistic/Poisson/ridge/AR(1) regression kernels.
- Added classical tests, one-way ANOVA and batched testing/correlation helpers.
- Added k-NN, graph/permutation algorithms and spatial median.
- Added directional-statistics and naive-Bayes modules.
- Added explicit 444-export coverage map; advanced routines not yet translated are marked pending instead of being silently omitted.
