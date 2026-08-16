# Porting notes

## Representation

`polynomial_t` stores an allocatable `real(dp)` coefficient vector in the same
increasing-power convention as the R package. Construction removes trailing
coefficients that are exactly zero by default. Numerical algorithms accept
explicit tolerances where rounding noise must be removed.

## Roots

The R package obtains roots from eigenvalues of a companion matrix. This port
uses a self-contained Aberth-Ehrlich simultaneous root iteration, avoiding a
BLAS/LAPACK dependency. Roots are sorted by real part and then imaginary part.
Exact zero roots are removed before iteration and restored afterward. Multiple
roots are intrinsically ill-conditioned; their returned values may form a
small cluster around the repeated root rather than identical bit patterns.

## Polynomial division and GCD

Long division is implemented directly. Remainders smaller than a scale-aware
floating-point tolerance are removed. This is necessary for practical GCD/LCM
calculations with floating-point coefficients. Exact operator equality remains
coefficient-wise.

## Orthogonal polynomials

The recurrence

```text
p_(j+1)(x) = (x-alpha_j) p_j(x) - beta_j p_(j-1)(x)
```

is constructed directly from values on the supplied abscissae, including
repeated abscissae as weights. With normalization enabled, evaluations have an
identity cross-product, matching `stats::poly` and upstream `poly.orth` up to
floating-point rounding and sign conventions. The leading coefficient remains
positive.

## R-specific behavior

- Plotting methods were omitted.
- S3 dispatch and runtime class checks were replaced by static derived types.
- R closures from `as.function.polynomial` were replaced by type-bound Horner
  evaluation.
- Polylist indexing, concatenation, repetition, and uniqueness map naturally to
  ordinary Fortran arrays and array assignment.
- Warning conditions are represented by `poly_status_t` codes/messages.
- Construction from complex roots is not exposed because upstream polynomials
  are real and upstream coercion discards imaginary coefficient parts.
