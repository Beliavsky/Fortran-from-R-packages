# API

## Main module

```fortran
use soma
```

The public floating-point kind is `dp = kind(1.0d0)`.

## Objective callback

```fortran
function objective(x) result(value)
    use soma, only : dp
    real(dp), intent(in) :: x(:)
    real(dp) :: value
end function objective
```

Additional problem data can be supplied through module state or a derived-type
wrapper in the caller, replacing R's `...` arguments.

## Bounds

```fortran
bnds = bounds(lower, upper)
```

`make_bounds` is an equivalent explicit name.

## Options

### All To One

```fortran
options = all2one(
    population_size=10,
    n_migrations=20,
    path_length=3.0_dp,
    step_length=0.11_dp,
    perturbation_chance=0.1_dp,
    min_absolute_sep=0.0_dp,
    min_relative_sep=1.0e-3_dp)
```

### T3A

```fortran
options = t3a(
    population_size=30,
    n_migrations=20,
    n_steps=45,
    migrant_pool_size=10,
    leader_pool_size=10,
    n_migrants=4,
    min_absolute_sep=0.0_dp,
    min_relative_sep=1.0e-3_dp)
```

### Pareto

```fortran
options = pareto(
    population_size=100,
    n_migrations=20,
    n_steps=10,
    perturbation_frequency=1.0_dp,
    step_frequency=1.0_dp,
    min_absolute_sep=0.0_dp,
    min_relative_sep=1.0e-3_dp)
```

## Optimization

```fortran
call soma_optimize(objective, bnds, result)
call soma_optimize(objective, bnds, result, options)
call soma_optimize(objective, bnds, result, options, init)
```

`init` has shape `(n_parameters, population_size)`, matching the R package's
column-oriented population representation.

The objective is a required procedure dummy with an explicit abstract
interface. The translated library does not call host-associated optional
procedure dummies from nested procedures, avoiding the gfortran
`-Werror=implicit-interface` portability problem encountered in some earlier
translations.

## Reproducible seeding

```fortran
call soma_set_seed(1234)
```

This seeds the Fortran intrinsic random-number generator deterministically for
a given compiler/runtime. It is not compatible with R's RNG stream.
