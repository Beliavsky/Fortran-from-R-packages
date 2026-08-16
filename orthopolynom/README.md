# orthopolynom-fortran

A modern Fortran 2018 translation of the computational code in the R package
`orthopolynom` 1.0-6.1, using the supplied `polynom-fortran` translation as its
polynomial arithmetic dependency.

## Implemented families

For each applicable family the library provides recurrence coefficients,
polynomials, normalized polynomials, norm-squared inner products, and weight
functions:

- Chebyshev C and S on `(-2,2)`
- Chebyshev T and U on `(-1,1)`
- shifted Chebyshev T and U on `(0,1)`
- Legendre and shifted Legendre
- spherical (Legendre alias used by the R package)
- Hermite H and probabilists' Hermite He
- generalized Hermite
- Laguerre and generalized Laguerre
- Gegenbauer / ultraspherical
- Jacobi P
- Jacobi G on `(0,1)`

## Generic computational API

The port also includes the non-family algorithms exported by `orthopolynom`:

- `pochhammer`, `lpochhammer`
- `orthogonal_polynomials` and `orthonormal_polynomials`
- `monic_polynomial_recurrences` and `monic_polynomials`
- `jacobi_matrices`
- `polynomial_coefficients`
- `polynomial_derivatives`
- `polynomial_integrals`
- `polynomial_orders`
- `polynomial_powers`
- `polynomial_roots`
- `polynomial_values`
- `polynomial_functions`
- `scale_x`

The R package's list of closures from `polynomial.functions` is represented by
`polynomial_function_list_t`; each item has an `evaluate` type-bound procedure.

## Example

```fortran
program example
  use orthopolynom
  implicit none
  type(polylist_t) :: p

  p = legendre_polynomials(4)
  print '(a)', p%item(5)%to_string()
  print '(f12.8)', p%item(5)%evaluate(0.25_dp)
end program example
```

## FPM

The supplied `polynom-fortran-v0.1.0` is vendored under `vendor/` and declared
as a local FPM dependency, so the release is self-contained.

```text
fpm test
fpm run --example demo_orthopolynom
```

## Validation

The tests cover:

- known coefficients for the classical and shifted families;
- recurrence and weight values;
- monic recurrences and Jacobi matrices;
- roots and polynomial-list helper functions;
- change-of-basis coefficients from `polynomial_powers`;
- finite and infinite-bound `scale_x` behavior;
- independent reference values for Gegenbauer, Jacobi P, Jacobi G,
  generalized Laguerre, and generalized Hermite;
- consistency of normalized and unnormalized parameterized families.

The release is tested with GNU Fortran using strict Fortran 2018 warnings and
runtime checking, and separately with optimization enabled.

## License

The translated `orthopolynom` code is GPL-2.0-or-later, matching the upstream
package's `GPL (>= 2)` declaration. The supplied `polynom-fortran` dependency
is GPL-2.0-only. See `LICENSES.md` and `NOTICE.md` for details. The complete
upstream `orthopolynom` source and archive are retained under `upstream/`.
