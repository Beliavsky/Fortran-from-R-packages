# API

All public facilities can be imported from module `mcga`.

## Kinds

- `dp` -- `real64` kind used throughout the package.
- `i32` -- integer kind used for byte values 0 through 255.

## Main optimizers

### `mcga_optimize`

```fortran
call mcga_optimize(popsize, chsize, minval, maxval, objective, result, &
                   crossprob, mutateprob, elitism, maxiter, seed)
```

`objective` has the explicit interface

```fortran
function objective(x) result(f)
  real(dp), intent(in) :: x(:)
  real(dp) :: f
end function objective
```

The returned `type(mcga_result)` contains

- `population(popsize,chsize)` sorted in increasing objective value;
- `costs(popsize)`;
- `generations`;
- type-bound `best()` returning the first chromosome.

`minval` and `maxval` initialize the population only. They are deliberately not
hard parameter constraints, matching upstream `mcga()`.

### `multi_mcga_optimize`

```fortran
call multi_mcga_optimize(popsize, chsize, numfunc, minval, maxval, &
                         objective, result, crossprob, mutateprob, &
                         elitism, maxiter, seed)
```

The callback is

```fortran
subroutine objective(x, f)
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: f(:)
end subroutine objective
```

`type(multi_mcga_result)` contains

- `population(popsize,chsize)`;
- `costs(popsize,numfunc)`;
- `ranks(popsize)` using the upstream rank-score definition;
- `generations`.

## Byte conversion

- `max_double()`
- `size_of_double()`
- `size_of_int()`
- `size_of_long()`
- `double_to_bytes(x)`
- `double_vector_to_bytes(x)`
- `bytes_to_double(bytes)`
- `byte_vector_to_doubles(bytes)`

Byte arrays use `integer(i32)` values in `[0,255]`.

## Byte genetic operators

- `one_point_crossover`
- `one_point_crossover_doubles`
- `two_point_crossover`
- `two_point_crossover_doubles`
- `uniform_crossover`
- `uniform_crossover_doubles`
- `byte_code_mutation`
- `byte_code_mutation_doubles`
- `byte_code_mutation_doubles_random`
- `ensure_bounds`

## GA-style operator library

These routines translate the computational operators from `R/oplibrary.R`:

- `byte_mutation`
- `byte_mutation_dynamic`
- `byte_mutation_random`
- `byte_mutation_random_dynamic`
- `byte_crossover`
- `byte_crossover_1p`
- `byte_crossover_2p`
- `sbx_crossover`
- `flat_crossover`
- `arithmetic_crossover`
- `blx_crossover`
- `linear_crossover`
- `unfair_average_crossover`

`linear_crossover` takes a `procedure(fitness_fn)` callback because upstream
selects the two offspring with highest GA fitness.

## Multi-objective rank score

`calculate_rank_scores(costs,ranks)` exposes the exact scoring loop used by
`multi_mcga.c`. It should not be confused with Pareto-front index or NSGA-II
nondomination rank.
