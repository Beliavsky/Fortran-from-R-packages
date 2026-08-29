# Porting notes

## Design

The R package accepts formulas and many fitted-model classes. The Fortran port uses explicit matrices and vectors, which makes the statistical kernels reusable without reproducing an R object system.

The least-squares layer uses LAPACK `DGELS`; covariance matrices use Cholesky inversion of cross-product matrices. Current high-level regression diagnostics therefore expect full-column-rank design matrices, matching the identifiable-model case used by the tests.

## Durbin-Watson

`src/pan.f` is the only compiled numerical source in upstream lmtest. It implements the amended Applied Statistics Algorithm AS 153 / AS R52. `lmtest_pan.f90` is a direct free-format module translation, retaining the algorithm and comments/provenance through the copy in `upstream/src/pan.f`.

For the exact Durbin-Watson path, the Fortran implementation forms the residual projection and a symmetric projected Durbin-Watson matrix. Its positive eigenvalues are the nonzero eigenvalues needed by AS 153. The large-sample fallback follows the mean/variance approximation in upstream `dwtest.R`.

## Distribution functions

The R implementation delegates tail probabilities to R's distribution functions. This port supplies the required normal, incomplete-beta, and incomplete-gamma calculations internally. Survival probabilities are evaluated directly where possible to avoid catastrophic cancellation in small p-values.

## R-specific features omitted

- Formula/model-frame construction and factor expansion.
- S3 methods and arbitrary fitted-model introspection.
- Automatic formula/model updating in `lrtest` and `waldtest`.
- `zoo` time indexes in `grangertest`; callers provide aligned numeric series.
- Printing/ANOVA-table decoration.
- `hmctest(plot=TRUE)` graphics.
- Bundled data sets and vignette rendering.

The upstream R sources are retained under `upstream/R` so these boundaries can be audited directly.
