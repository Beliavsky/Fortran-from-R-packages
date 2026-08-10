# Translation coverage

## Translated

- R `nnls()` numerical behavior.
- R `nnnpls()` numerical behavior.
- Lawson-Hanson active/passive-set iteration.
- Mixed-sign column transformation used by `nnnpls`.
- `MODE` status semantics and the default `3*N` iteration limit.
- `NSETP`, passive/bound index sets, residual norm, residuals, fitted values,
  and deviance.

## Modernization choices

The original package calls fixed-form Fortran 77 routines that update a QR
factorization in place using `H12` Householder transformations and `G1`
rotations. The modern implementation retains the Lawson-Hanson active-set
logic but recomputes the passive-set least-squares problem with a modern
Householder QR solve whenever the passive set changes. This removes GOTOs,
implicit interfaces, transformed input work arrays, and old storage-stride
APIs while preserving the mathematical algorithm.

Near-dependent candidate columns are rejected using a scale-aware QR rank
threshold. This serves the same purpose as the original `DIFF` independence
test, although the exact branch can differ on pathological nearly singular
inputs.

## Omitted R-only code

- S3 classes and print methods.
- R input checks for NA/NaN/NULL/infinite values. Fortran callers should pass
  finite numeric arrays.
- R accessor generics (`coef`, `fitted`, `residuals`, `deviance`). Their data
  are fields of `nnls_result` directly.

## Numerical fidelity

A separate validation program linked both the modern implementation and the
original fixed-form Lawson-Hanson routines. Across 60 deterministic random
matrix problems, testing both NNLS and NNNPLS (120 solves total), the largest
absolute difference in any coefficient or residual norm was about 1.1e-14.
