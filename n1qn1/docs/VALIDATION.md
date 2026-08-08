# Validation

## Supplied Rosenbrock regression

For the package's three-dimensional Rosenbrock-like objective starting from
`[1.02, 1.02, 1.02]`, the Fortran translation reproduces:

- final parameters `[1, 1, 1]`
- final objective `1`
- 40 objective and gradient evaluations
- the supplied reference Hessian within `2e-9` maximum absolute error

With the documented prior Hessian, it reproduces 29 evaluations and the
reference final Hessian within `2e-9`.

Reusing the dense Hessian from the first fit through `initial_hessian` produces
33 evaluations, matching the supplied package regression.

## Additional tests

- diagonal positive-definite quadratic
- single-variable quadratic
- exact internal factor restart
- evaluation-limit termination
- polymorphic user data
- progress callback cancellation
- strict explicit-interface compilation
- runtime bounds and floating-point checking
