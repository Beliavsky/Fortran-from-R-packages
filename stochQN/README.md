# stochQN-fortran

A modern Fortran 2018 translation of the computational code in the R package
`stochQN` 0.1.2-1 by David Cortes.

The library implements stochastic limited-memory quasi-Newton methods for
smooth optimization:

- online L-BFGS (`oLBFGS`)
- stochastic quasi-Newton (`SQN`)
- adaptive quasi-Newton (`adaQN`)

The translation preserves the original reverse-communication design. An
optimizer returns a typed request describing the calculation that the caller
must provide next: a regular gradient, a same-batch gradient, a large-batch
gradient, a Hessian-vector product, or a validation objective value.

Higher-level callback runners and a native stochastic logistic-regression model
are also included.

## Build with FPM

```text
fpm build
fpm test
fpm run --example quadratic_sqn
fpm run --example stochastic_logistic
```

The package has no external numerical dependencies. Dot products, matrix-vector
products, norms, and elementwise operations use standard Fortran intrinsics.

## Main modules

- `stochqn_kinds`: double-precision kind
- `stochqn_core`: reverse-communication optimizer types and status constants
- `stochqn_guided`: callback-driven complete optimization runs
- `stochqn_logistic`: logistic loss, gradient, Hessian-vector product, prediction,
  and an online typed model

## Low-level example

```fortran
use stochqn_kinds, only : dp
use stochqn_core

type(olbfgs_t) :: optimizer
type(stochqn_request_t) :: request
real(dp) :: x(2), gradient(2)

x = [2.0_dp, -1.0_dp]
gradient = 0.0_dp
call optimizer%initialize(2, mem_size=10)
call optimizer%advance(0.1_dp, x, gradient, request)

do while (optimizer%get_iteration() < 100 .or. request%task /= task_calc_grad)
   select case (request%task)
   case (task_calc_grad, task_calc_grad_same_batch)
      gradient = request%x
   case default
      error stop 'Unexpected request.'
   end select
   call optimizer%advance(0.1_dp, x, gradient, request)
end do
```

## Numerical and interface differences

The numerical state machines and correction-pair formulas follow the supplied C
kernel. R formulas, S3 objects, data frames, factor handling, printing, and R's
reference-semantics wrappers are not reproduced. They are replaced by typed
Fortran objects and explicit array interfaces.

See `docs/API_MAP.md`, `docs/PORTING_NOTES.md`, and `docs/VALIDATION.md` for
full details.

## License

BSD 2-Clause, matching the original package. The complete original package
source is retained under `original/stochQN-master` for provenance.
