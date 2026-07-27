# Porting notes

## Language mapping

SharpeR's S3 classes are represented by explicit derived types. Numeric vectors
and matrices are accepted directly, optional arguments replace R defaults, and
status fields report recoverable numerical failures. Procedures are exposed
through the umbrella module `sharper`.

R's column-major matrix semantics align naturally with Fortran. One-dimensional
Fortran arrays replace R vectors, while observations occupy rows and assets
occupy columns in return matrices.

## Numerical distributions

The project is self-contained. It does not call R, SciPy, Boost, GSL, LAPACK,
or another statistics package.

- Noncentral t probabilities and densities use adaptive numerical integration
  of the normal/chi-square representation.
- Noncentral F probabilities and densities use a centered Poisson mixture of
  incomplete-beta terms.
- Noncentral chi-square probabilities use a Poisson mixture of incomplete-gamma
  terms.
- Quantiles use safeguarded bracketing and bisection.
- Random generators use Fortran's intrinsic uniform generator plus native
  normal, gamma, chi-square, t, and F transformations.

These routines prioritize portability and accuracy over the throughput of a
specialized external distribution library.

## Conditional maximum test

The R implementation delegates polyhedral conditional-normal calculations to
`epsiwal`. The Fortran interface instead accepts the selected contrast,
polyhedral constraint matrix, bounds, covariance matrix, and observed vector
explicitly, then computes the one-dimensional truncated-normal probability
natively. This preserves the computational test while removing the R package
adapter.

## Unified covariance estimators

The Gaussian branch of `sm_vcov` is implemented from the analytic covariance
identities for sample first and second moments. It is mathematically equivalent
to the normal-model asymptotic calculation but is not a line-for-line port of
R's `matrixcalc` construction. The empirical branch follows the stacked-moment
sample-covariance definition. `ism_vcov` propagates uncertainty with a central
finite-difference Jacobian. After the delta-method matrix products, the routine
averages each pair of opposite-triangle entries and copies that common value to
both positions. This enforces the exact symmetry required of a covariance matrix
and avoids compiler-dependent roundoff differences observed with GNU Fortran on
Windows.

## Documented corrections

Two indexing/scaling issues in the supplied R source are not reproduced:

1. In the hedged Markowitz helper, the supplied source overwrites the mean with
   its transformed value and later maps `transpose(H)` against that transformed
   mean. The dimensions and intended optimization imply that the solved
   transformed portfolio weights should be mapped back instead. The Fortran
   implementation returns `matmul(transpose(H), work_w)`.
2. The supplied unpaired-test path can annualize values twice when the source
   objects and requested output use nonunit observations per epoch. The Fortran
   implementation first converts each estimate and standard error to native
   epoch units, forms the contrast once, and annualizes the result once.

Both choices are covered by deterministic tests and are called out here rather
than silently presented as literal source reproduction.

## Approximation boundaries

- Maximum-Sharpe confidence handling uses the translated chi-square, Follman,
  and Bonferroni approximations rather than R object-specific root wrappers.
- Very extreme noncentral parameters may require more computation than mature
  platform libraries because all tail calculations are performed natively.
- Singular covariance matrices return a nonzero status instead of invoking an
  R pseudoinverse fallback.

## Excluded surfaces

Printing, plotting, serialized package data, formula evaluation, `lm` extraction,
`xts`/`zoo`/`timeSeries` conversion, and vignette generation are excluded. They
are environment integrations rather than computational Sharpe-ratio methods.
