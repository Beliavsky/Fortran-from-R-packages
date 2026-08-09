# Translation coverage

## Translated

The complete computational body of upstream `R/dykstra.R` is represented:

- argument/default semantics;
- diagonal versus general-matrix paths;
- `factorized=FALSE` and `factorized=TRUE` transformations;
- positive-semidefinite spectral regularization;
- unconstrained starting solution;
- epsilon rescaling by the unconstrained solution magnitude;
- per-constraint Dykstra correction vectors;
- cyclic equality/inequality projections;
- convergence and iteration-limit behavior;
- retransformation to the original coordinates;
- quadratic objective evaluation when the original `D` is available.

R's `eigen(..., symmetric=TRUE)` is replaced by a self-contained cyclic Jacobi
symmetric eigensolver.  Eigenvalues/eigenvectors are sorted in descending
order to match the assumptions made by the R code.

## R-only code omitted

- `print.dykstra` formatting;
- R list/class creation;
- R coercion and S3 infrastructure.

The equivalent information is represented by `type(dykstra_result)`.

## Source quirks and deliberate handling

### Diagonal `unconstrained` field

For a diagonal `D`, Dykstra 1.0-0 computes the reported `unconstrained` field
as

```text
dvec / Rinv^2
```

where `Rinv = 1/sqrt(diag(D))`.  Algebraically, the true unconstrained
minimizer would instead use `Rinv^2*dvec`.  The Fortran translation preserves
the upstream expression in `result%unconstrained` for source compatibility.
The optimized `result%solution` is retransformed correctly and is not affected
by this reporting quirk.

### Singular PSD matrices

The upstream code shifts small eigenvalues, but the all-zero diagonal case can
still leave a zero shift and lead to division by zero.  The Fortran port uses
`max(1, spectral_scale)` when forming the regularization threshold and ensures
a strictly positive shift.  This is a defensive extension of the documented
upstream intent (support positive-semidefinite `D`).

### Zero constraint columns

A violated zero constraint column would divide by zero in the R implementation.
The Fortran routine returns a negative status instead.

### Exact equality test

The R routine tests equality constraints with exact `Ai == b[i]`.  To remain
warning-clean under strict gfortran without introducing a numerical tolerance,
the Fortran code expresses exact equality as simultaneous `Ai <= b(i)` and
`Ai >= b(i)`.
