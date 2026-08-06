# minqa-fortran

A modern Fortran 2018 port of the computational code in the R package
`minqa` 1.2.8.

The package provides Powell's derivative-free optimization algorithms based
on quadratic interpolation:

- `bobyqa`: bound-constrained optimization
- `newuoa`: unconstrained optimization with a flexible interpolation set
- `uobyqa`: unconstrained optimization with a full quadratic interpolation set

The original fixed-form Fortran kernels were converted to free-form source,
placed in a module with explicit interfaces, and connected to a native
Fortran procedure callback. The R, Rcpp, S3, and dynamic-library layers are
not required.

## Build with FPM

```text
fpm build
fpm test
fpm run --example minqa_example
```

## Minimal example

```fortran
program example
   use minqa_module, only : dp, minqa_control_t, minqa_result_t, bobyqa
   implicit none

   real(dp) :: x(2)
   type(minqa_control_t) :: control
   type(minqa_result_t) :: result

   x = [-1.2_dp, 1.0_dp]
   control%rhobeg = 0.25_dp
   control%rhoend = 1.0e-7_dp

   call bobyqa(rosenbrock, x, result, &
      [-2.0_dp, -1.0_dp], [2.0_dp, 3.0_dp], control)

   write(*, '("x =",*(1x,f12.8))') result%x
   write(*, '("f =",es14.6)') result%fval

contains

   function rosenbrock(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f

      f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
   end function rosenbrock

end program example
```

## Public API

The module `minqa_module` exports:

- `dp`
- `objective_function`
- `minqa_control_t`
- `minqa_result_t`
- `bobyqa`
- `newuoa`
- `uobyqa`
- `minqa_status_message`

`x` is updated in place and is also copied into `result%x`.

## Status codes

| Status | Meaning |
|---:|---|
| 0 | normal exit |
| 1 | maximum objective evaluations exceeded |
| 2 | invalid interpolation-point count |
| 3 | trust-region step failed |
| 4 | bound interval too small |
| 5 | excessive cancellation in an update denominator |
| 6 | invalid public-API input |
| 7 | unknown core status |

`result%raw_status` retains the original Powell/minqa status code.

## Notes

The native callback bridge uses one module-level active procedure pointer
during an optimization call. Therefore, calls are not reentrant or
thread-safe within one process. Independent processes are unaffected.

See `docs/API_MAP.md`, `docs/PORTING_NOTES.md`, and
`docs/VALIDATION.md` for details.

## License

GPL-2.0-only, matching the source R package. Original source files are
retained under `original/minqa-master` for attribution and traceability.
