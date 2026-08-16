# Porting notes

## Scope

This port translates the computational R code exported by `orthopolynom`
1.0-6.1. The supplied package contains no plotting implementation to omit.
R data frames and lists are represented by Fortran derived types and allocatable
arrays.

The supplied `polynom-fortran-v0.1.0` translation is used rather than
reimplementing polynomial arithmetic. It is included unchanged under
`vendor/polynom-fortran` and referenced as a local FPM dependency.

## Representation

A recurrence table is represented by `recurrence_t`, containing arrays `c`,
`d`, `e`, and `f`. The package recurrence convention is

```text
c_j p_(j+1)(x) = (d_j + e_j x) p_j(x) - f_j p_(j-1)(x).
```

`monic_recurrence_t` contains the corresponding `a` and `b` arrays. Variable
length vectors and matrices are represented by list derived types in
`orthopolynom_types`.

## `polynomial.functions`

R can return lists of closures. Fortran has no direct equivalent with the same
lightweight ownership semantics, so the port returns a
`polynomial_function_list_t`. Each item owns a `polynomial_t` and provides
scalar/vector `evaluate` methods.

## `polynomial.roots`

The R code forms symmetric Jacobi matrices and calls R's `eigen`. The port
still exposes `jacobi_matrices`, but computes `polynomial_roots` from the
corresponding monic polynomials using the supplied `polynom-fortran` complex
root solver. Orthogonal-polynomial roots are real; the returned list contains
the real parts, sorted in ascending order.

## Numerical changes

Several formulas using factorials are evaluated through `log_gamma` where this
is algebraically equivalent, which delays overflow.

`scaleX` in the upstream R source contains an `is.finite(u) || v == Inf`
condition that overwrites the intended `(-Inf, Inf)` result. The Fortran port
implements the evident four intended cases:

- `(-Inf, Inf)`: unchanged;
- `(-Inf, v)`: upper-end translation;
- `(u, Inf)`: lower-end translation;
- finite `(u,v)`: affine rescaling.

The upstream `lpochhammer(z, 0)` returns `1` even though the mathematical
logarithm of `(z)_0 = 1` is zero. The Fortran port intentionally preserves the
upstream return value for compatibility.

The upstream `ultraspherical.weight` evaluates `log(1-x^2)` without checking
support, unlike `gegenbauer.weight`, even though ultraspherical is otherwise
an alias for Gegenbauer. The Fortran port returns zero outside `(-1,1)`,
consistent with the rest of the package's weight functions.

The special `jacobi.g` `p=0` recurrence code writes a second row even when
`n=0`. The Fortran implementation guards that row so degree-zero requests are
valid without out-of-bounds access.
