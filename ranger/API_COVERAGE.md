# API and algorithm coverage

This document maps the computational surface of upstream `ranger` 0.18.0 to the native Fortran API. R formula parsing, S3 methods, printing, data-frame conversion, sparse-R-object adapters, and plotting/presentation code are intentionally not translated.

## Forest engines

| Upstream area | Fortran status | Native entry points / notes |
|---|---|---|
| Classification forest | Translated | `fit_ranger_classification`, `predict_ranger_classification`; Gini, Hellinger, ExtraTrees; OOB votes/error |
| Probability forest | Translated | `fit_ranger_probability`, `predict_ranger_probability`; class-probability terminal estimates and OOB probability loss |
| Regression forest | Translated | `fit_ranger_regression`, `predict_ranger_regression`; variance, ExtraTrees, Maxstat, Beta and Poisson split criteria |
| Survival forest | Translated | `fit_ranger_survival`, `predict_ranger_survival`; log-rank, ExtraTrees, AUC/AUC-ignore-ties and Maxstat paths; Nelson–Aalen CHF; `time.interest` explicit grids are sorted/deduplicated and scalar counts select approximately equally spaced observed event times |
| Quantile regression forest | Translated computational core | `options%quantreg` prepares random terminal-node response values for ordinary and OOB prediction, matching ranger quantile preparation; `predict_ranger_quantiles` takes quantiles across treewise node values |
| Terminal node IDs / predict-all | Translated | optional `terminal_nodes`, `per_tree`, `per_tree_probability`, `per_tree_chf` outputs |

## Sampling, predictor handling, and tree controls

- sampling with and without replacement;
- scalar sample fraction and class-wise sample fractions for classification/probability;
- nonnegative case weights, including zero-weight holdout cases;
- user-supplied in-bag count matrices;
- `mtry`, `min.node.size`, `min.bucket`, and `max.depth`; classification/probability also accept class-specific node/bucket vectors, including zero entries exactly as ranger does;
- numeric and integer-coded categorical predictors;
- unordered-factor ordered and partition searches;
- learned missing-value left/right routing;
- split-selection weights and always-split variables;
- ExtraTrees random split candidates;
- regularization factors, including depth-dependent regularization;
- class weights;
- in-bag retention;
- `options%oob_error` independently disables cumulative OOB prediction/error work while leaving permutation-importance computation available.

## Importance and inference

| Upstream area | Status | Notes |
|---|---|---|
| Impurity importance | Translated | split-decrease accumulation |
| Permutation importance | Translated | classification, probability, regression, survival |
| Casewise permutation importance | Translated | stored as predictor × observation matrix |
| Scaled permutation importance | Translated | standard-error scaling from treewise importance changes |
| Corrected impurity importance | Translated | A forest-wide shadow predictor block is built with one shuffled row permutation shared across all predictors, as in ranger's `Data::permuteSampleIDs`; shadow variables participate in candidate weighting/always-split logic and their signed impurity decreases are accumulated back to the corresponding original variables. |
| Janitza importance p-values | Translated | `importance_pvalues_janitza` mirrors non-positive importance values as upstream R code does |
| Altmann p-values | Computational test translated | `importance_pvalues_altmann` computes p-values from an observed vector and a caller-generated null-importance matrix. Formula/data-frame response permutation and repeated R object construction are omitted. |
| Infinitesimal jackknife | Translated with optimizer adaptation | Raw IJ, Monte Carlo bias correction, the without-replacement finite-sample correction, ranger's `n <= 20` no-calibration behavior, and the empirical-Bayes `gfit`/`gbayes` convolution/interpolation workflow are translated. The two-parameter minimization uses a native Newton/line-search routine instead of R's `nlm`, so optimizer trajectories are not claimed bit-for-bit identical. |
| Maxstat significance correction | Translated | Native average ranks/logrank scores, cutpoint accounting, Lau92 and Lau94 corrections, Benjamini-Hochberg adjustment, and the survival one-cut unadjusted-normal special case mirror ranger's regression/survival paths. |

## Higher-level computational helpers

- hierarchical shrinkage for regression forests: translated;
- hierarchical shrinkage for probability forests: translated;
- case-specific random-forest proximity weights: translated as `case_specific_weights`; callers can use these weights with the normal fit routines to perform the repeated case-specific fitting step;
- holdout random forests: represented directly by `options%holdout=.true.` plus zero/nonzero case weights; the R wrapper that randomly constructs two folds is orchestration rather than a distinct forest algorithm;
- tree sizes and variable-use counts: translated.

## Adaptations relative to C++ implementation details

The statistical algorithms are translated, but C++-specific infrastructure is not carried over literally:

- Rcpp/RcppEigen and R external-pointer/object serialization;
- memory modes (`double`, `float`, compressed char) and SNP two-bit storage;
- sparse `dgCMatrix` input adapters;
- C++ thread scheduling/OpenMP execution framework;
- specialized small-Q/large-Q split-search caches and data sorting caches;
- binary forest serialization files.

The Fortran implementation uses direct candidate evaluation and dense native arrays. This can differ in speed and tie-breaking from the highly optimized C++ core while retaining the statistical criteria and forest semantics. Regression and probability OOB outputs use NaN for observations with zero OOB coverage, matching ranger; survival leaves uncovered CHF rows at zero. The native classification OOB class label is integer-valued, so zero is used as the no-coverage sentinel and `oob_count` is the authoritative coverage mask.

## R-only/interface code omitted

- formulas, terms, data-frame/factor discovery, `gwaa.data`, and R `Matrix` adapters;
- S3 printing/coercion/extraction methods (`predictions`, `importance`, `timepoints` are direct fields/native outputs instead);
- `treeInfo()` data-frame presentation and `deforest()` R object trimming;
- deprecated `getTerminalNodeIDs()` wrapper (terminal IDs are native prediction outputs);
- R warnings/messages/progress callbacks and seed integration with R's RNG;
- C++ thread scheduling and the `num.threads` control; `node_stats` is retained in `ranger_options`, but native tree objects always keep the node statistics used by the translated post-fit helpers instead of suppressing those fields for memory savings;
- R object serialization and prediction-object construction;
- arbitrary R `what=` functions for quantile prediction (native quantile probabilities are supported directly).
