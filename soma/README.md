# soma-fortran

Modern Fortran translation of the computational code in the R package
`soma` 1.2.0 by Jon Clayden.

The package implements the Self-Organising Migrating Algorithm (SOMA) with
all three strategies exported by the R package:

- All To One
- Team To Team Adaptive (T3A)
- Pareto

The translation is standalone and uses only standard Fortran. There is no R,
Rcpp, BLAS or LAPACK dependency.

## Build

With FPM:

```text
fpm build
fpm test
```

Two examples are included:

```text
fpm run --example rastrigin
fpm run --example strategies
```

A strict GNU Fortran validation script is available as
`scripts/test_gfortran.sh` and `scripts/test_gfortran.bat`.

## Basic use

```fortran
program example
    use soma, only : dp, soma_bounds, soma_result, bounds, soma_optimize, soma_set_seed
    implicit none

    type(soma_bounds) :: bnds
    type(soma_result) :: result

    interface
        function sphere(x) result(value)
            use soma, only : dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function sphere
    end interface

    bnds = bounds([-5.0_dp, -5.0_dp], [5.0_dp, 5.0_dp])
    call soma_set_seed(1234)
    call soma_optimize(sphere, bnds, result)

    print *, result%population(:,result%leader)
    print *, result%cost(result%leader)
end program example

function sphere(x) result(value)
    use soma, only : dp
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sum(x*x)
end function sphere
```

The option constructors intentionally retain the R package names:
`all2one()`, `t3a()` and `pareto()`.

## Result fields

`type(soma_result)` contains:

- `leader`: 1-based index of the final best individual
- `population`: final parameter vectors, one individual per column
- `cost`: final cost for each individual
- `history`: best cost at initialization and after each counted migration
- `evaluations`: evaluation-count history using the same convention as the R package
- `migrations`: number of counted migrations
- `status` and `message`: Fortran validation/status information

See `API.md` and `TRANSLATION_COVERAGE.md` for details.

## License

GPL-2.0-only, matching the upstream package. The original source package is
included under `original/soma-master/` for provenance.
