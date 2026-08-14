# API mapping

## Upstream package structure

The supplied mlr 2.19.3 tree exports 432 R names and contains 150
`RLearner_*` adapter files: 75 classification, 54 regression, 9 survival,
10 clustering and 2 multilabel adapters. Most adapters call algorithms from
Suggests/Imports packages and contain no implementation of those algorithms.
Their filenames are recorded in `upstream/EXTERNAL_LEARNER_ADAPTERS.txt`.

## Direct computational mappings

| mlr concept/function family | Fortran counterpart |
|---|---|
| regression measures | `mlr_metrics::measure_*` |
| classification measures | `mlr_metrics::measure_*` |
| `calculateConfusionMatrix` | `confusion_matrix` |
| `.632` aggregation | `aggregate_b632` |
| holdout / crossval / repeated CV | `make_holdout`, `make_kfold`, `make_repeated_kfold` |
| subsampling | `make_subsample` |
| bootstrap OOB | `make_bootstrap_oob` |
| `imputeMean`, `imputeMedian`, constants | `impute_mean`, `impute_median`, `impute_constant` |
| scaling | `fit_standardizer`, `apply_standardizer`, `standardize` |
| dummy features | `one_hot_encode` |
| down/oversampling | `downsample_classes`, `oversample_classes` |
| native SMOTE C kernel | `nearest_neighbors`, `smote_generate` |
| `classif.featureless` | `fit_featureless_classifier`, `predict_featureless_classifier` |
| `regr.featureless` | `fit_featureless_regression`, `predict_featureless_regression` |
| `regr.lm` | `fit_linear_regression`, `predict_linear_regression` |
| `classif.logreg` | `fit_logistic_regression`, logistic prediction routines |
| `cluster.kmeans` | `fit_kmeans`, `predict_kmeans` |
| k-NN-style learner | `fit_knn_*`, `predict_knn_*` |
| survival Cox learner | `fit_cox_learner`, `predict_cox_risk` |
| `cindex` | `measure_cindex` |
| Kaplan-Meier computations | `fit_kaplan_meier` |
| resample/evaluate | `resample_regression`, `resample_classification` |
| grid tuning | `grid_search` |
| random tuning | `random_search` |
| exhaustive feature selection | `feature_select_exhaustive` |
| sequential feature selection | `feature_select_forward` |
| random feature selection | `feature_select_random` |
| binary `tuneThreshold` | `tune_binary_threshold` |

## Dependency-bound upstream APIs

The `RLearner_*` files that merely register/call randomForest, ranger,
glmnet, xgboost, e1071, kernlab, nnet, mboost, gbm, H2O, RWeka, C50,
DiceKriging and many other packages are catalogued but not translated as mlr
algorithms. They should be connected to independent Fortran ports of those
packages as those become available.

Similarly, tuning backends whose optimization algorithms live in external
packages (`cmaes`, `GenSA`, `mlrMBO`, `irace`, NSGA-II backends) are represented
in v0.1.0 by generic grid/random callback search rather than mislabeled
reimplementations.

## R-framework-only omissions

Not applicable to a Fortran numerical library: S3 classes and printing,
formula/model-frame parsing, data.frame/data.table operations, package
registry/discovery, XML configuration, plotting, cache management,
parallelMap/batchtools orchestration, R error/dump objects, and vignette/UI
helpers.

## Future work

Reasonable v0.2 candidates include permutation feature importance, partial
dependence, calibration curves, Friedman/Nemenyi benchmark statistics,
multiclass threshold tuning, multilabel measures, `.632+`, additional filter
scores, bagging/stacking abstractions, and adapters to already translated
Fortran learner packages.
