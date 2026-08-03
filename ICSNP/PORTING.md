# Porting notes

## General mapping

The R package API was translated to array-oriented Fortran procedures. Dots in
R names become underscores. Optional R arguments become Fortran optional
arguments, and list/`htest` results become derived types.

## Linear algebra

The project contains self-contained implementations of pivoted dense linear
solves, matrix inversion and determinants, Jacobi symmetric eigendecomposition,
symmetric square roots, inverse square roots, covariance matrices, and
Mahalanobis distances.

The routines target small and medium multivariate-statistics problems. They do
not attempt to replace optimized BLAS/LAPACK for very large matrices.

## Robust estimators

`tyler_shape` uses the standard fixed-point covariance-form Tyler iteration,
normalizing each iterate to determinant one. The original R package sometimes
iterates on the inverse shape; the two forms have the same normalized fixed
point.

`duembgen_shape` applies Tyler estimation to all nonzero pairwise differences.
Weighted variants use products of observation weights.

The spatial median follows the modified Weiszfeld update used by the original
R code, including the coincident-observation correction.

## External R dependencies

The original package uses `ICS::ics` and `ICS::ics.components` inside
`ind.ictest`. This translation uses covariance whitening followed by a
FOBI-style fourth-moment eigendecomposition. It provides the same required
invariant-coordinate role without an external dependency, but component order
and signs can differ from the R `ICS` defaults.

The original package uses `mvtnorm::rmvnorm` for simulation-based calibration.
This translation includes a portable Box-Muller normal generator.

## Interfaces omitted

The following are R infrastructure rather than separate numerical algorithms:

- formula methods;
- S3 generic dispatch;
- `ics` object methods;
- construction and printing of R `htest` objects;
- `na.action` dispatch;
- package datasets and documentation-browser behavior.

No plotting routines are present in the original package.

## Numerical distributions

Chi-square and F survival probabilities, normal quantiles, and chi-square
quantiles are implemented internally using incomplete gamma/beta continued
fractions, a rational normal-quantile approximation, and safeguarded bisection.
