# subplex-fortran

Modern Fortran translation of the computational code in the R package
`subplex` 1.9.

The package implements Tom Rowan's Subplex algorithm for unconstrained
minimization. Subplex repeatedly partitions the coordinates according to the
magnitudes of recent changes and applies Nelder-Mead simplex minimization on
each low-dimensional subspace.

## Build

```sh
fpm build
fpm test
```

A strict GNU Fortran validation script is available as
`scripts/test_gfortran.sh` and `scripts/test_gfortran.bat`.

## Basic use

```fortran
use subplex, only: dp, subplex_result, subplex_minimize

type(subplex_result) :: result
real(dp) :: x0(2)

x0 = [11.0_dp, -33.0_dp]
call subplex_minimize(objective, x0, result)
```

See `API.md` and `example/` for complete examples.

## License and provenance

The upstream R package declares GPL-3. The full GPL text supplied by upstream
is copied to `LICENSE`, and the unmodified uploaded package is retained under
`original/subplex-master/`.
