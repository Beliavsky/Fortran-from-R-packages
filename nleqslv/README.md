# nleqslv-fortran

A standalone modern Fortran/FPM port of the computational core of the R package
`nleqslv` 3.3.7.

The package solves square systems of nonlinear equations using Newton or
Broyden iterations with the globalization strategies provided by upstream
`nleqslv`:

- cubic line search
- quadratic line search
- geometric line search
- double-dogleg trust region
- Powell single-dogleg trust region
- More-Hebden / Levenberg-Marquardt hook step
- pure Newton/Broyden step (`none`)

The public interface is Fortran 2018 and uses typed option/result derived types
and procedure callbacks. The upstream fixed-form solver kernels were converted
to free-form `.f90` source and detached from R/C. LAPACK and BLAS remain the
only external numerical-library requirements. Version 0.1.1 also supplies
explicit interfaces for all BLAS/LAPACK and internal solver calls, so strict
FPM builds with `-Werror=implicit-interface` are supported.

## Build

With FPM and system LAPACK/BLAS:

```text
fpm build
fpm test
fpm run --example
```

The manifest links `lapack` and `blas`.

## Basic use

```fortran
program demo
   use nleqslv_fortran
   implicit none
   type(nleq_options) :: opt
   type(nleq_result) :: sol
   real(dp) :: x0(2)

   x0 = [2.0_dp, 0.5_dp]
   opt = nleq_options()
   opt%method = NLEQ_BROYDEN
   opt%global = NLEQ_DBLDOG

   call solve_nleqslv(x0, equations, sol, opt)
   print *, sol%x
   print *, trim(sol%message)

contains
   subroutine equations(x, f)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
      f(1) = x(1)**2 + x(2)**2 - 2.0_dp
      f(2) = exp(x(1)-1.0_dp) + x(2)**3 - 2.0_dp
   end subroutine equations
end program demo
```

A user Jacobian can be supplied as the final argument to `solve_nleqslv`.
For numerical banded Jacobians, set `dsub` and `dsuper` in `nleq_options`.

## Public API

- `solve_nleqslv` -- counterpart of R `nleqslv`
- `search_zeros` -- counterpart of R `searchZeros`
- `test_nleqslv` -- computational counterpart of R `testnslv`
- `termination_message`
- `nleq_options`, `nleq_result`, `search_zeros_result`, `nleq_test_result`
- method/global/scaling integer constants

See `API_MAP.md` and `PORTING_NOTES.md` for details.

## Dependencies

System BLAS and LAPACK are required. No R runtime or C shim is required.

## License

GPL-2.0-or-later, following upstream `nleqslv`. See `LICENSES.md`, `NOTICE.md`,
and the retained source tree under `upstream/`.
