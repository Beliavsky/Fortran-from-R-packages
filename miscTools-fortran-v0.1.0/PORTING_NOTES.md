# Porting notes

## Scope

This port targets the computational content of `miscTools` 0.6-30. Plotting,
console-reporting/digest code, arbitrary R attributes, and S3 fitted-model
dispatch are not reproduced.

## Semidefiniteness

The upstream package offers two algorithms:

- `method="det"`: all principal minors;
- `method="eigen"`: symmetric eigenvalues.

Both are available. The determinant path enumerates principal submatrices and
uses pivoted Gaussian elimination. The eigen path uses a standalone Jacobi
symmetric-eigensolver, so the package has no LAPACK dependency.

The default follows upstream behavior: determinant checking below dimension 13
and eigenvalue checking for dimension 13 or larger.

Matrices are accepted if symmetric to `1000*tol` and are explicitly
symmetrized before testing, matching the R source.

## Medians and R NA values

Fortran IEEE NaN is used as the analogue of R `NA` in the median routines.
With `na_rm=.false.`, a NaN input yields NaN; with `na_rm=.true.`, NaNs are
removed. Matrix, rank-3, and rank-4 `colMedians` variants are supplied, covering
the upstream tests and common array use.

## Coefficient tables

The upstream `coefTable` calls `pt` for two-sided Student-t p-values. The
Fortran port evaluates the equivalent regularized-incomplete-beta expression
locally, avoiding an external statistics dependency.

When no degrees of freedom are supplied, the fourth column uses `-1` as a
Fortran sentinel instead of R `NA`.

## Matrix/vector ordering

`sym_matrix`, `vecli`, `vecli2m`, `triang`, and `veclipos` preserve R's
column-major triangular ordering. In particular:

`vecli([[11,12,13],[12,22,23],[13,23,33]])`
returns `[11,12,13,22,23,33]`.

## R object interfaces

`stdEr` in R can discover `std`, `vcov`, and coefficient names through generic
model objects. Fortran instead exposes the numerical kernel directly:
`std_er(cov,se)`.

Similarly, `nObs`/`nParam` are represented by numerical matrix/vector helpers;
there is no attempt to emulate arbitrary R S3 model structures.
