# Translation coverage

## Directly translated computational behavior

- `marqLevAlg` / `mla` optimization loop.
- `deriva` numerical first derivatives and information matrix.
- `deriva_grad` numerical information matrix from an analytical gradient.
- `searpas`, `valfpa`, `func1`, and `ghg` logic, absorbed into module
  procedures.
- Adaptive diagonal inflation and positive-definiteness retries.
- Full and partial-Hessian RDM convergence.
- `loglikLMM` and `gradLMM`.

## Modernized native kernels

The upstream `dchole.f90`, `dmfsd.f90`, and `dsinv.f90` use packed storage,
GOTO-heavy Fortran, and legacy implementation details.  The translation keeps
packed upper-triangle interoperability but performs the required SPD solve and
inverse with modern dense Cholesky routines.

For `marqLevAlg`, this is algorithmically equivalent at the decision points:
`dchole` is used only to determine whether the inflated information matrix is
positive definite and, once it is, solve for the step; `dsinv` is used to
invert a positive-definite information matrix for RDM/covariance calculations.
Floating-point roundoff can differ from the legacy kernels.

## Deliberately omitted R infrastructure

- `foreach`, `doParallel`, and `parallel::makeCluster` orchestration.
- S3 `print`/`summary` methods.
- R call/model objects, `...`, package/export registration, and `.Fortran`
  registration glue.
- Console/file iteration reporting.
- Example `.rda` serialization.

The Fortran finite-difference implementation is serial.  Parallel evaluation in
the R package changes execution scheduling, not the mathematical derivative
formulas.

## Source quirk preserved

Although the R documentation calls `hess` a Hessian callback, the source copies
its output directly into the matrix subsequently treated as the information
matrix.  The Fortran overload does the same; it does not silently negate or
rescale a user-supplied matrix.
