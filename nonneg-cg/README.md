# nonneg-cg-fortran

Modern Fortran translation of the computational kernel in the R package
`nonneg.cg` 0.1.6-1.

The package minimizes a differentiable objective subject to componentwise
nonnegativity using the modified Polak-Ribiere-Polyak conjugate-gradient method
implemented by David Cortes from Li (2013).

## Build

```text
fpm build
fpm test
fpm run --example rosenbrock
```

The package is self-contained and does not require BLAS.

## Minimal use

```fortran
use nonneg_cg
real(dp) :: x(2)
type(nonneg_cg_result_t) :: result

x = [0.0_dp, 2.0_dp]
call minimize_nonneg_cg(x, objective, gradient, result)
```

Callbacks may optionally receive polymorphic user data. A monitor callback can
cancel a run after a completed iteration.

## Evaluation counts

`result%nfeval` preserves the original C package's reported counter. The C code
increments this counter only after rejected line-search trials, so it is not the
true number of objective calls. `result%objective_calls` reports the actual
number of objective evaluations.

## License

BSD-2-Clause, preserving the original copyright and disclaimer.
