# ppcor-fortran

A self-contained modern Fortran translation of the computational code in the
R package **ppcor 1.1**, which calculates partial and semi-partial (part)
correlations and their significance tests.

The port provides Pearson, Spearman, and Kendall methods, pairwise test
wrappers, tied-rank handling, Kendall tau-b, Student-t and normal p-values, and
a Moore-Penrose pseudoinverse fallback for rank-deficient association matrices.
It has no external numerical-library dependency.

## Build with FPM

```text
fpm build
fpm test
fpm run
```

Run an individual example with, for example:

```text
fpm run --example example_partial_matrix
```

## Basic use

```fortran
use ppcor, only : dp, ppcor_result, pcor, ppcor_pearson

real(dp) :: x(100,4)
type(ppcor_result) :: result

! Fill x with observations in rows and variables in columns.
call pcor(x, result, ppcor_pearson)
if (result%status == 0) then
   print *, result%estimate
   print *, result%p_value
end if
```

The public aliases `pcor` and `spcor` correspond to the R functions. The more
descriptive names `partial_correlation` and `semi_partial_correlation` are also
available. `pcor_test` and `spcor_test` accept either a vector or matrix of
conditioning variables.

## Matrix orientation

Input is `observations x variables`, matching the ordinary R matrix layout.
Result matrices are `variables x variables`.

Semi-partial correlations are directional and therefore generally asymmetric.
For `estimate(i,j)`, the other variables are removed from variable `i`, matching
the orientation produced by the upstream matrix formula.

## Numerical behavior

For Pearson data the upstream code starts from a covariance matrix, while this
port starts from a correlation matrix. Partial and semi-partial correlations are
invariant to positive rescaling, so the resulting coefficients are equivalent.
Spearman uses average ranks for ties. Kendall uses tau-b.

When the association matrix is numerically rank deficient, a symmetric
Moore-Penrose pseudoinverse is formed by Jacobi eigendecomposition. The result
record reports `used_pseudoinverse` and the estimated rank.

## Scope

The complete computational API of the upstream package is represented. R data
frames, list objects, vector recycling, warnings, and package documentation
machinery are replaced by typed Fortran arrays, result records, and status codes.
Missing/non-finite values are rejected explicitly rather than propagated through
R's `NA` rules.

See [API_MAP.md](API_MAP.md) and [PORTING_NOTES.md](PORTING_NOTES.md).

## License

GPL-2.0-only, matching the upstream package. The complete upstream source
snapshot is retained under `upstream/ppcor-master`.
