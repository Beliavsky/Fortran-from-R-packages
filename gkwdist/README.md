# gkwdist-fortran

Modern Fortran translation of the computational code in R package `gkwdist`
1.1.4, with an FPM project layout.

## Implemented families

- generalized Kumaraswamy (GKw)
- beta-Kumaraswamy (BKw)
- Kumaraswamy-Kumaraswamy (KKw)
- exponentiated Kumaraswamy (EKw)
- McDonald
- Kumaraswamy
- the package's beta parameterization, Beta(gamma, delta+1)

For every family the library provides density, CDF, quantile, RNG, negative
log-likelihood, analytical gradient, and analytical Hessian routines. It also
provides `gkwgetstartvalues`, the package's multi-start method-of-moments
initialization algorithm.

## Example

```fortran
program demo
   use gkwdist
   implicit none
   real(dp) :: p

   p = pgkw(0.4_dp, 1.7_dp, 2.4_dp, 1.3_dp, 0.8_dp, 1.2_dp)
   print *, p
   print *, qgkw(p, 1.7_dp, 2.4_dp, 1.3_dp, 0.8_dp, 1.2_dp)
end program demo
```

The `d/p/q` routines are elemental and therefore work naturally with
conformable arrays and scalar parameters.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The supplied `numDeriv-fortran` translation is vendored as a dependency and is
used by `test_derivatives` to independently check the analytical derivatives.
The production gkwdist library does not call numDeriv.

## Tests

- `test_distributions`: closed-form and independent reference values, tails,
  quantile inversion, nested-family identities
- `test_derivatives`: analytical gradients/Hessians versus numDeriv
- `test_rng`: support and beta mean check
- `test_startvalues`: method-of-moments starting-value smoke tests
- `test_api`: exercises all 49 family d/p/q/r/ll/gr/hs exports

See `API_MAP.md` and `PORTING_NOTES.md` for translation details.
