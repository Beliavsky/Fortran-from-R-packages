# maxLik-fortran

Modern Fortran 2018 translation of the computational core of the R package
`maxLik` 1.5-2.2.

The library provides a common callback-based interface for maximum-likelihood
and general nonlinear maximization. It includes deterministic, derivative-free,
and stochastic optimizers; numerical derivative tools; fixed parameters;
linear constraints; and likelihood-based covariance diagnostics.

## Build with FPM

```text
fpm build
fpm test
fpm run
```

The public module is:

```fortran
use maxlik
```

No external numerical library is required.

## Minimal example

```fortran
program example
  use maxlik, only: dp, maxlik_problem, maxlik_result, &
    initialize_problem, max_lik
  implicit none

  type(maxlik_problem) :: problem
  type(maxlik_result) :: result

  call initialize_problem(problem, 2, objective)
  problem%gradient => gradient
  call max_lik(problem, [-3.0_dp, 4.0_dp], result, 'bfgs')

  print *, result%estimate
  print *, result%maximum

contains

  subroutine objective(x, value, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    value = -0.5_dp * ((x(1) - 1.0_dp)**2 + (x(2) + 2.0_dp)**2)
    status = 0
  end subroutine objective

  subroutine gradient(x, g, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(out) :: status
    g = [-(x(1) - 1.0_dp), -(x(2) + 2.0_dp)]
    status = 0
  end subroutine gradient

end program example
```

A callback returns `status = 0` on success. Analytic gradients, Hessians, and
observation-level scores are optional. Missing derivatives are computed
numerically. Observation-level scores are required by BHHH, SGA, and Adam.

## Implemented scope

- Unified `max_lik` dispatcher and common result/control types.
- Newton-Raphson with step halving and Marquardt-style Hessian correction.
- BFGS and BFGSR maximization.
- BHHH/Fisher-scoring updates from observation-level scores.
- Nonlinear conjugate-gradient and Nelder-Mead maximization.
- Simulated annealing with deterministic seeded random-number generation.
- Stochastic-gradient ascent with momentum and Adam.
- Central or forward numerical gradients, numerical Hessians, and derivative
  comparison reports.
- Numerical Jacobians for vector-valued callbacks.
- Active/fixed parameters and box bounds.
- Linear equality and inequality constraints in upstream form
  `A * theta + b = 0` and `A * theta + b >= 0`, implemented by an outer
  quadratic-penalty sequence.
- Hessian covariance, observation-score robust covariance, standard errors,
  normal confidence intervals, AIC, and condition-number diagnostics.
- Stored objective values and parameter paths when requested in the control.
- Active-parameter packing and expansion utilities corresponding to the role
  of upstream subset/fixed-parameter helpers.

## Interface differences from R

R functions accept closures, `...`, S3/S4 objects, named vectors, and objective
values carrying gradient/Hessian attributes. Fortran uses explicit procedure
pointers and arrays. Supply callbacks through a `maxlik_problem`; results are
returned in a `maxlik_result` derived type.

The translated constraints use a quadratic-penalty sequence rather than R's
`constrOptim` barrier machinery. The reported `constraint_violation` should be
checked after constrained fitting. Fixed-parameter covariance rows and columns
are zero, matching the upstream convention.

R printing, `summary`, `tidy`, `glance`, formula/model-frame handling, object
attributes, package reexports, and vignette generation are not computational
algorithms and are not translated. The package contains no plotting code.

See `API.md`, `PORTING.md`, `TRANSLATION_COVERAGE.md`, and `TESTING.md`.

## License

The upstream package is licensed under GPL version 2 or any later version. This
translation preserves that license and is distributed under
`GPL-2.0-or-later`. Full GPL version 2 and version 3 texts are included. The
complete upstream package snapshot is retained under `original/` for attribution
and traceability.
