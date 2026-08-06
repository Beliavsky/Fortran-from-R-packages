# polynom-fortran

A modern Fortran 2018 translation of the computational code in the R package
`polynom` 1.4-1.

The library represents a real univariate polynomial by coefficients in
increasing powers:

```text
coef(1) + coef(2) x + ... + coef(n) x^(n-1)
```

## Features

- Polynomial construction with trailing-zero removal
- `+`, `-`, `*`, `/`, `**`, `==`, and `/=` operator overloads
- Polynomial quotient and remainder
- Scalar and vector Horner evaluation
- Derivatives, antiderivatives, and definite integrals
- Lagrange interpolation and construction from real roots
- Origin shifts, monic normalization, GCD, and LCM
- Complex roots using a self-contained Aberth iteration
- Zeros, stationary points, and inflexion-point summaries
- Orthogonal or orthonormal polynomials on weighted abscissae
- Polylist sums, products, derivatives, integrals, GCD, and LCM
- Coefficient rounding, floor, ceiling, truncation, and significant digits
- Compatibility names corresponding closely to the R API

Plotting and R's S3 object, expression, and printing infrastructure are not
included. `polynomial_t%to_string()` supplies reusable text formatting, while
ordinary Fortran arrays replace R list indexing, concatenation, repetition,
and uniqueness operations.

## Build

With GNU Make and GNU Fortran:

```sh
make checked
make optimized
```

With FPM:

```sh
fpm test
fpm run --example demo_polynom
```

The checked build uses strict Fortran 2018 conformance, warnings as errors,
bounds/runtime checking, and backtraces.

## Minimal example

```fortran
program example
  use polynom
  implicit none
  type(polynomial_t) :: p, dpdx

  p = poly_from_roots([1.0_dp, 2.0_dp, 3.0_dp])
  dpdx = derivative(p)

  print '(a)', p%to_string()
  print '(a,es16.8)', 'p(4) = ', p%evaluate(4.0_dp)
  print '(a)', dpdx%to_string()
end program example
```

## License

GPL-2.0-only, matching the upstream package. The complete upstream source and
original ZIP are retained under `upstream/`.
