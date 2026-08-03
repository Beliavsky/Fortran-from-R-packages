# Porting notes

## Interface mapping

R matrices map directly to rank-two `real(dp)` arrays. R lists and classed return
objects are replaced by Fortran derived result types. Optional R arguments map to
optional Fortran arguments or fields of `oracle_control`. R function names are
converted to lower-case underscore names.

The translation does not include formula, data-frame, plotting, or package-data
interfaces. These are not part of the numerical library.

## Linear algebra

RcppArmadillo operations are replaced by self-contained dense Gaussian
elimination, Jacobi symmetric eigendecomposition, thin SVD, Cholesky, and
pseudoinverse routines. These are suitable for small and medium research
problems but do not match optimized BLAS/LAPACK throughput.

Near-singular systems use a symmetric pseudoinverse where mathematically
appropriate and return explicit status values instead of terminating R.

## HAC covariance

The implementation uses the Newey-West lag rule
`floor(4*(n/100)^(2/9))`, Bartlett weights, symmetric positive/negative lag
contributions, and optional VAR(1) prewhitening with covariance reversion.
This also avoids the scalar pass-by-value and one-sided accumulation hazards in
some upstream native-code paths.

## Identification tests

`iterative_kleibergen_paap_2006_beta_rank_test` is an asymptotic
singular-value implementation based on scaled factor loadings. It preserves the
iterative rank-testing API but is not the full covariance-weighted
Kleibergen-Paap statistic used by the R package. The result message identifies
this adaptation. The Chen-Fang bootstrap is implemented around the translated
scaled-loading/SVD machinery and therefore inherits this difference when the KP
rank estimate initializes the bootstrap.

## Oracle TFRP

The closed-form adaptive soft threshold is retained. GCV, k-fold CV, rolling
validation, the one-standard-error rule, relaxed refitting, and the three
upstream weighting choices are supported. Folds are deterministic and
sequential. No parallel backend is used.

## FGX testing

The `glmnet` dependency is replaced by a deterministic coordinate-descent Lasso
with an internally generated lambda path and modulo fold assignment. Therefore,
selected controls can differ from R `glmnet` because standardization, lambda
construction, convergence paths, and fold assignment are not identical.
The final covariance is a conventional HAC sandwich based on selected controls
and new-factor residuals.

## Numerical and reproducibility notes

- Random bootstrap operations use an internal xorshift generator.
- Supply explicit 64-bit seeds for reproducible Chen-Fang results.
- All calculations are double precision.
- Input must be finite; missing-value imputation is outside this package.
- Standard errors can differ slightly across compilers because dense
  eigendecomposition and pseudoinverse thresholds are portable implementations.
