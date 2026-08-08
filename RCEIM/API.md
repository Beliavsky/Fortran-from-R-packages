# API

## `module rceim`

### `type(rceim_options)`

Fields:

- `n_total` - population size, default 1000
- `n_elite` - elite count; `<=0` uses `n_total/4`
- `n_super` - elite members copied unchanged, default 1
- `alpha` - CE location smoothing coefficient, default 1
- `epsilon` - positive fractional convergence tolerance, default 0.1
- `q` - dynamic sigma-smoothing exponent, default 2
- `max_iter` - maximum iteration counter, default 50
- `wait_gen` - generations without best-score change; `<=0` uses `max_iter`
- `chaos_gen` - generations before chaos restart; `<=0` uses `max_iter`
- `minimize` - `.true.` for minimization, `.false.` for maximization
- `verbose` - print per-iteration best/mean/sd
- `seed` - nonnegative value seeds the Fortran RNG

### `type(rceim_result)`

- `x(:)` - best parameter vector
- `value` - objective value in its natural sign
- `score` - internal minimized score (`value` when minimizing, `-value` when maximizing)
- `converged` - convergence flag
- `iterations` - final iteration counter, matching upstream loop semantics
- `criterion` - stopping explanation
- `elite_x(:,:)` - final elite parameter matrix
- `elite_score(:)` - final elite internal scores
- `history_best(:)`, `history_mean(:)`, `history_sd(:)` - elite score histories
- `history_n` - number of populated history entries

### `ceim_optimize(fn, lower, upper, result [, options])`

Runs the translated CE-inspired optimizer. `fn` has the explicit abstract interface `rceim_objective`.

### `enforce_domain(params, lower, upper)`

Clips every parameter to its nearest boundary, matching `enforceDomainOnParameters`.

### `sort_population_by_score(params, scores)`

Stable ascending sort of population rows by score, which is the computational use of upstream `sortDataFrame` inside `ceimOpt`.

## `module rceim_benchmarks`

- `test_fun_optimization(x)`
- `test_fun_optimization_2d(x)`

These are direct translations of the two exported R benchmark functions.

## `module rceim_random`

- `rceim_set_seed(seed_value)`
- `fill_uniform(x, lo, hi)`
- `fill_normal(x, mu, sigma)`
- `random_normal_scalar()`
