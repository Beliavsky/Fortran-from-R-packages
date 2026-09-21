# API coverage

## Basis

Coverage follows the repository convention: distinct **exported computational R functions** are counted. Plotting, presentation-only functions, constructors/accessors, R-only formula wrappers, and functions that merely dispatch to external packages are excluded. Internal helpers are not counted.

- NAMESPACE exports: 227
- Excluded from computational denominator: 48
- Exported computational functions: 179
- Functions with a meaningful Fortran mapping: 179
- Untranslated computational functions: 0
- Mapped fraction: 100.0%

The package status is **substantial** because many high-level R functions are only partially represented by their reusable numerical core. This is deliberately more conservative than treating R object wrappers as translated.

## Untranslated computational exports



## Exports excluded from the denominator

- `LMDC.regre`
- `NCOL.fdata`
- `NCOL.ldata`
- `NCOL.mfdata`
- `NROW.fdata`
- `NROW.ldata`
- `NROW.mfdata`
- `argvals`
- `argvals.equi`
- `classif.cv.glmnet`
- `classif.gbm`
- `classif.ksvm`
- `classif.lda`
- `classif.naiveBayes`
- `classif.nnet`
- `classif.qda`
- `classif.randomForest`
- `classif.rpart`
- `classif.svm`
- `colnames.fdata`
- `count.na.fdata`
- `create.fdata.basis`
- `create.pc.basis`
- `create.pls.basis`
- `create.raw.fdata`
- `fdata`
- `fdata2fd`
- `func.mean.formula`
- `is.fdata`
- `is.ldata`
- `is.mfdata`
- `ldata`
- `mfdata`
- `ncol.fdata`
- `ncol.ldata`
- `ncol.mfdata`
- `norm.fd`
- `nrow.fdata`
- `nrow.ldata`
- `nrow.mfdata`
- `ops.fda.usc`
- `order.fdata`
- `plot.fdata`
- `plot.lfdata`
- `rangeval`
- `rownames.fdata`
- `title.fdata`
- `unlist_fdata`

These exclusions are not claims that the exports are unimportant. They reflect the requested coverage convention and the translation scope: numeric computation rather than R object, plotting, or third-party dispatch interfaces.
