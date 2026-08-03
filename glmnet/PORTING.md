# Porting notes

## What was translated

The computational model families and exported numerical utilities from
`glmnet` 5.0 were translated to modern Fortran. The project retains the
upstream R and C++ source under `original/` for provenance.

## Main numerical adaptations

### Common path engine

The upstream package uses a large family of specialized templated C++ kernels,
strong rules, active sets, and sparse iterators. This port uses:

- cyclic coordinate descent for Gaussian and IRLS weighted least squares;
- block coordinate descent for multiresponse Gaussian;
- proximal gradient with backtracking for multinomial and Cox models;
- warm starts across decreasing lambda values.

The objectives and penalties are retained, but iteration counts and the last
few digits of coefficients need not match the C++ engine.

### Sparse input

`glmnet_sparse_csc` preserves an explicit CSC interface. The current fitting
wrapper materializes the matrix as dense storage before fitting. It is useful
for interoperability, not for reproducing upstream sparse scalability.

### Cox regression

The Cox implementation supports right-censored and `(start, stop]` data,
strata, observation weights, offsets, Breslow ties, and Efron ties. Risk sets
are evaluated directly, giving approximately quadratic work in the number of
observations for general data. The upstream engine is substantially faster.

The returned Cox `dev_ratio` is improvement in average negative partial
log-likelihood relative to the zero-coefficient model. It is not guaranteed to
be bit-identical to the R package's partial-deviance normalization.

### General GLM families

R family objects cannot cross a Fortran boundary directly. The typed
`fit_custom_family_path` callback provides the corresponding IRLS extension:
it returns a working response, IRLS weights, and deviance for a supplied
linear predictor.

### Relaxed fits

Active sets are refitted with a tiny ridge (`lambda=1e-10`, `alpha=0`) instead
of an exactly singular unpenalized solve. This is stable when active columns
are collinear or when an active set is wider than the effective sample rank.

### Cross-validation

Cross-validation is sequential and deterministic. R `foreach` parallelism is
not reproduced. Out-of-fold predictions are stored directly in
`glmnet_cv_result%predictions`.

### Data preparation

`makeX` and `prepareX` in R understand data frames and factor contrasts. The
Fortran forms operate on numeric matrices, replace nonfinite values, remove
constant columns, and concatenate train/test matrices. Categorical encoding
must be supplied explicitly by the caller.

## Omitted R infrastructure

The compiled library omits:

- S3 printing and plotting
- formula/model-frame processing
- survival `survfit` objects
- R progress bars and `foreach` backends
- R Matrix class dispatch
- bundled `.rda` datasets and vignette build products

These omissions do not remove a numerical model family from the Fortran API.
