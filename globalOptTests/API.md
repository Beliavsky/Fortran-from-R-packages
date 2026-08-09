# API

Module: `global_opt_tests`

## Kinds and catalogue

```fortran
integer, parameter :: dp
integer, parameter :: n_benchmarks = 50
character(len=16), parameter :: benchmark_names(n_benchmarks)
```

## Generic benchmark dispatcher

```fortran
real(dp) function go_test(x, fn_name, check_dim, status)
```

`check_dim` defaults to `.true.`.  Status values are:

- `0`: success
- `1`: wrong parameter dimension
- `2`: unknown benchmark name

On an error, the returned value is quiet NaN.

## Metadata

```fortran
subroutine get_default_bounds(fn_name, lower, upper, status)
integer function get_problem_dimension(fn_name, status)
real(dp) function get_global_opt(fn_name, status)
```

Bounds are allocated by `get_default_bounds`.

## Direct objective procedures

All 50 upstream objective names are available as pure scalar functions with
signature

```fortran
pure real(dp) function name(x)
    real(dp), intent(in) :: x(:)
```

The public functions are:

`ackleys`, `aluffipentini`, `beckerlago`, `bohachevsky1`,
`bohachevsky2`, `branin`, `camel3`, `camel6`, `cosmix2`, `cosmix4`,
`dekkersaarts`, `easom`, `emichalewicz`, `expo`, `goldprice`,
`griewank`, `gulf`, `hartman3`, `hartman6`, `hosaki`, `kowalik`, `lm1`,
`lm2n10`, `lm2n5`, `mccormic`, `meyerroth`, `mielecantrell`,
`modlangerman`, `modrosenbrock`, `multigauss`, `neumaier2`, `neumaier3`,
`paviani`, `periodic`, `powellq`, `pricetransistor`, `rastrigin`,
`rosenbrock`, `salomon`, `schaffer1`, `schaffer2`, `schubert`,
`schwefel`, `shekel10`, `shekel5`, `shekel7`, `shekelfox5`, `wood`,
`zeldasine10`, and `zeldasine20`.
