# Rsolnp-fortran

A modern Fortran translation of the computational core of the R package
`Rsolnp` 2.0.1. The library solves smooth nonlinear optimization problems
with box constraints, equality constraints, and two-sided nonlinear
inequality constraints by an augmented-Lagrangian method with projected BFGS
inner iterations.

The project is self-contained and uses no external numerical libraries.

## Implemented functionality

- `solnp` and `csolnp` single-start constrained optimization
- Optional analytic objective gradients and constraint Jacobians
- Central finite-difference derivatives when analytic callbacks are absent
- Equality, two-sided inequality, and parameter-bound constraints
- Feasibility restoration, adaptive penalty updates, and BFGS curvature updates
- KKT diagnostics and regularized least-squares multiplier estimation
- Deterministic low-discrepancy starting-point generation
- Sequential multistart optimization through `csolnp_ms` and `gosolnp`
- Standard-form problem conversion
- Complete 77-entry benchmark registry and 18 executable benchmark definitions

The R printing methods, parallel backend, vignettes, and package-specific list/S3
infrastructure are not part of the Fortran library.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example basic_constrained
fpm run --example multistart_example
fpm run --example standard_form_example
```

The package uses standard Fortran 2018 and has no FPM dependencies.

## Minimal use

```fortran
use rsolnp, only : dp, solnp_problem, solnp_result, csolnp

type(solnp_problem) :: problem
type(solnp_result) :: result

problem%n = 2
problem%fn => objective
allocate(problem%start(2), problem%lower(2), problem%upper(2))
problem%start = [-1.2_dp, 1.0_dp]
problem%lower = -5.0_dp
problem%upper = 5.0_dp

call csolnp(problem, result)
```

Callback signatures and complete type descriptions are in [API.md](API.md).

## Numerical and API differences from R

Fortran callbacks and derived types replace R functions, lists, S3 objects, and
`...` arguments. `solnp` and `csolnp` share the translated solver; `csolnp`
uses supplied analytic derivatives and otherwise falls back to finite
differences. Multistart evaluation is deterministic and sequential. The
returned Hessian is the final BFGS approximation for the augmented objective.
See [PORTING.md](PORTING.md) for details.

## License

The upstream package declares GPL-2. This translation is distributed under
**GPL-2.0-only**. Original metadata and computational sources are retained under
`original/`.
