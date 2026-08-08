# API

## Precision

`dp = kind(1.0d0)` is exported by module `oor`.

## POO

```fortran
call poo(f, horizon, noise_level, result [, rhomax, nu])
```

Arguments:

- `f`: scalar objective callback; POO maximizes it on `[0,1]`.
- `horizon`: nominal function-evaluation budget.
- `noise_level`: bound used to form `alpha = log(horizon)*noise_level**2`.
- `rhomax`: number of equidistant rho values in `[0,1]`, default 20.
- `nu`: HOO metric scale, default 1.

`type(poo_result)` contains `par`, `value`, `best_rho`, evaluation/sample
counts, and the final `poo_tree`.

## StoSOO / SOO

```fortran
call stosoo(fn, lower, upper, nb_iter, result [, options])
```

`type(stosoo_options)` fields:

- `stochastic = .true.`: `.false.` selects deterministic SOO.
- `maximize = .false.`: the R API minimizes by default.
- `keep_tree = .false.`: retain the normalized search tree in the result.
- `verbose = 0`: retained for API compatibility; the library itself is quiet.
- `k_max = 0`: zero requests the upstream-derived default.
- `h_max = 0`: zero requests the upstream-derived default.
- `delta = -1`: nonpositive requests `1/sqrt(nb_iter)`.

`type(stosoo_result)` contains:

- `par`: best point in the original user bounds.
- `value`: objective value using the user's min/max convention.
- `evaluations`: number of objective evaluations actually performed.
- `xs(:,i)`, `ys(i)`: evaluation history in original coordinates/sign.
- `tree`: optional normalized tree when `keep_tree=.true.`.

## Benchmark functions

The module exports `guirland`, `sin1`, `difficult`, `difficult2`, and
`double_sine`.

At the singular source points, `difficult(0)` and `difficult2(0.5)` return a
quiet NaN rather than raising an R-language condition.
