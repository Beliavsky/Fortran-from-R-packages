# API

## Main module

```fortran
use trust
```

This re-exports the public API from the implementation modules.

## `trust_options`

Fields:

- `rinit`: initial trust-region radius. Default `1`.
- `rmax`: maximum trust-region radius. Default `100`.
- `iterlim`: maximum trust-region subproblems. Default `100`.
- `fterm`: termination tolerance for actual objective change.
- `mterm`: termination tolerance for predicted model change.
- `minimize`: `.true.` for minimization, `.false.` for maximization.
- `save_history`: retain the upstream `blather=TRUE` style iteration history.
- `parscale(:)`: optional positive finite parameter scaling vector.

## Objective interface

```fortran
subroutine trust_objective(x, value, gradient, hessian, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    real(dp), intent(out) :: gradient(:)
    real(dp), intent(out) :: hessian(:, :)
    integer, intent(out) :: status
end subroutine
```

For finite objective values, gradient and Hessian must also be finite.  An infeasible trial point may be signaled with `+Inf` during minimization or `-Inf` during maximization.  A nonzero `status` signals evaluation failure.

## `trust_optimize`

```fortran
call trust_optimize(objective, parinit, options, result)
```

`result` contains:

- `argument(:)`
- `value`
- `gradient(:)`
- `hessian(:,:)`
- `converged`
- `iterations`
- `status`
- `message`
- `history` when enabled

## Step types

Integer constants are provided for the exact upstream trust-region subproblem classifications:

- `trust_step_newton`
- `trust_step_easy_easy`
- `trust_step_hard_easy`
- `trust_step_hard_hard`

`trust_step_name(code)` returns the corresponding upstream-style text label.

## History

When `save_history=.true.`, `result%history` contains:

- `argument(:,i)` -- current iterate before subproblem `i`
- `argument_try(:,i)` -- proposed point
- `step_type(i)`
- `accepted(i)`
- `radius(i)`
- `step_norm(i)`
- `rho(i)`
- `value(i)`
- `value_try(i)`
- `predicted_difference(i)`

Only entries `1:history%n` are meaningful.
