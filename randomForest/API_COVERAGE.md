# API and computational coverage

This translation targets the numerical behavior owned by `randomForest`; it does not attempt to reproduce R objects, formulas, factors, plotting devices, S3 dispatch, or data-frame metadata.

## Translated computational areas

| Upstream area | Fortran API / status |
| --- | --- |
| classification forest engine (`randomForest.default`, `rf.c`, `classTree.c`, `rfsub.f`) | `fit_classification` |
| classification prediction | `predict_classification` |
| numeric CART Gini splits | translated |
| categorical CART Gini splits | exhaustive small-category and randomized large-category search translated |
| class weights and cutoffs | translated |
| weighted bootstrap / sampling without replacement | translated |
| stratified sampling | translated |
| OOB votes and class/OOB error curves | translated |
| mean decrease in accuracy and Gini impurity | translated |
| OOB/all-case proximities and in-bag counts | translated |
| regression engine (`regrf.c`, `regTree.c`) | `fit_regression` |
| regression prediction and OOB MSE | `predict_regression`, forest fields |
| numeric/categorical regression splits | translated |
| regression permutation/impurity importance | translated |
| optional linear bias correction | translated |
| unsupervised synthetic-class forest | `fit_unsupervised`, `make_synthetic_class` |
| `margin` | `classification_margin` |
| `outlier` | `outlier_scores` |
| numeric `na.roughfix` | `roughfix_numeric` |
| `classCenter` numerical center calculation | `class_centers` |
| `treesize` | generic `tree_sizes` |
| `varUsed` | generic `variable_usage` |
| `partialPlot` numerical averaging | `partial_dependence_classification`, `partial_dependence_regression` |
| `rfImpute` numerical iterations | `rf_impute_classification`, `rf_impute_regression` |
| `tuneRF` numerical OOB search | `tune_classification_mtry`, `tune_regression_mtry` |
| `rfcv` numerical feature-ranking CV | `rfcv_classification`, `rfcv_regression` |
| `MDSplot` numerical `cmdscale(1-proximity)` step | `mds_coordinates`; plotting omitted |

## Intentional omissions or interface differences

- `randomForest.formula`, factor/data-frame conversion, R `terms`, S3 methods, printing, plotting, `varImpPlot`, and `rfNews` are R-specific interfaces.
- `MDSplot` graphics are omitted; only the classical MDS coordinates/eigenvalues are exposed.
- `partialPlot` graphics are omitted; only the numerical partial-dependence curve is exposed.
- `combine`, `grow`, and `getTree` primarily manipulate or display R forest objects. Fortran forest components and individual `rf_tree` values are directly accessible as typed derived types.
- Test-set bookkeeping embedded in the R fit object is replaced by explicit prediction calls.
- `localImp` case-wise importance is not currently stored as a forest field; global class/overall permutation importance and its standard error are translated.
- `rfImpute` accepts a numeric `real(dp)` matrix with optional integer-coded categorical columns instead of a mixed-type R data frame.
- `rfcv` exposes deterministic Fortran fold generation and feature-count/error/prediction arrays; R names and list structure are omitted.

## Indexing and categorical representation

Fortran observations and variables use one-based indexing naturally. Categorical predictor levels are integer codes `1:ncat(j)` represented exactly in a `real(dp)` column. `ncat(j)=1` means numeric. Trees store categorical left-branch membership explicitly in logical masks rather than exposing the upstream packed double bit mask.
