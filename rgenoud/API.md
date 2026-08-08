# API

## `genoud_optimize`

```fortran
call genoud_optimize(fn, lower, upper, options, result, &
                     gradient, starting_values, hessian)
```

`fn(x)` returns a scalar objective. `lower` and `upper` are length-`nvars`
bounds. Optional `gradient(x,g)` supplies an analytical Euclidean gradient.
`starting_values` is an `nvars x nstart` matrix. Set `hessian=.true.` to return
a central-difference Hessian in `result%hessian`.

## `genoud_optimize_lexical`

```fortran
call genoud_optimize_lexical(fn, nfit, lower, upper, options, result, &
                             comparator, local_objective, local_gradient, &
                             starting_values, hessian)
```

`fn(x,f)` fills `f(1:nfit)`. Without a comparator, the criteria are compared
lexicographically, all in the minimization or maximization direction selected
by `options%maximize`. A custom `comparator(a,b)` can implement any strict
ordering.

As in upstream rgenoud, derivative refinement of a lexical problem needs a
scalar local objective. Pass it as `local_objective`; otherwise BFGS, gradient
checking, and operator 9 are disabled for that run.

## `genoud_options`

Important fields mirror the R interface:

- `pop_size = 1000`
- `max_generations = 100`
- `wait_generations = 10`
- `hard_generation_limit = .true.`
- `maximize = .false.`
- `memory_matrix = .true.`
- `solution_tolerance = 1e-3`
- `gradient_check = .true.`
- `use_bfgs = .true.`
- `boundary_enforcement = 0`
- `integer_parameters = .false.`
- `operator_weights(1:9) = 50,...,50,0`
- `p9_mix = -1` (random mixing weight)
- `bfgs_burnin = 0`
- `seed`
- `print_level`

Operator indices correspond to R's `P1` through `P9`:

1. cloning/survival
2. uniform mutation
3. boundary mutation
4. non-uniform mutation
5. polytope crossover
6. simple crossover
7. whole non-uniform mutation
8. heuristic crossover
9. local-minimum crossover

## `genoud_result`

- `par(:)` best parameters
- `fit(:)` best scalar or lexical fitness vector
- `gradient(:)` gradient at the returned solution when applicable
- `hessian(:,:)` optional numerical Hessian
- `generations`
- `peak_generation`
- `pop_size`
- `operators(1:9)` realized operator population counts
- `evaluations`
- `unique_evaluations`
- `converged`
- `status`

## Derivative/statistics helpers

The public module also exports `numerical_gradient`, `numerical_hessian`, and
`sample_moments`.
