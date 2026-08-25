# numDeriv-fortran

A modern Fortran translation of the computational code in the R package
`numDeriv` 2016.8-1.1. The library calculates numerical gradients,
Jacobians, Hessians, and Bates-Watts derivative matrices.

## Features

- Richardson extrapolation with configurable step sizes and iteration count
- Simple forward or backward first differences
- Per-coordinate one-sided Richardson derivatives
- Complex-step gradients and Jacobians
- Hybrid complex-step/Richardson Hessians
- Bates-Watts `D` matrices for scalar- or vector-valued functions
- Explicit status codes for invalid options, nonfinite values, and shape errors
- No external numerical dependencies

## Build

```text
fpm build
fpm test
fpm run
```

The package version in `fpm.toml` is `2016.8.1`, a valid semantic-version
representation of the upstream version `2016.8-1.1`.

A direct GNU Fortran validation script is also provided:

```text
./scripts/test_gfortran.sh
```

On Windows:

```text
scripts\test_gfortran.bat
```

## Minimal example

```fortran
module functions
   use numderiv, only : dp
   implicit none
contains
   function objective(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sum(sin(x))
   end function objective
end module functions

program example
   use functions, only : objective
   use numderiv, only : dp, grad
   implicit none
   real(dp) :: x(3), g(3)

   x = [0.0_dp, 0.5_dp, 1.0_dp]
   call grad(objective, x, g)
   print *, g
end program example
```

See `API.md`, `PORTING.md`, and the `example/` directory for the complete
callback interfaces and complex-step variants.

## Scope

All four exported computational routines from the R package are represented:
`grad`, `jacobian`, `hessian`, and `genD`. R S3 dispatch and R-specific object
printing are not applicable to Fortran and are omitted.

## License

GPL-2.0-or-later, following the upstream package's `GPL-2` license declaration.
The upstream computational sources and metadata are retained under `original/`.
