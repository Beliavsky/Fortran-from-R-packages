# API summary

All public BB routines are re-exported by:

```fortran
use bb
```

## Optimization

```fortran
fit = spg(par, fn [, control] [, gr])
fit = spg_box(par, fn, lower, upper [, control] [, gr])
fit = spg_projected(par, fn, project [, control] [, gr])
fit = spg_linear(par, fn, a, b, meq [, control] [, gr] [, lower] [, upper])
```

`spg_result` fields are `par`, `value`, `gradient`, `fn_reduction`, `iter`,
`feval`, `convergence`, `method`, `m`, and `message`.

The SPG convergence codes follow the R package convention:

- 0: success
- 1: maximum iterations
- 2: maximum function evaluations
- 3: function-evaluation failure
- 4: gradient-evaluation failure
- 5: projection/search-direction failure

The robust driver is:

```fortran
fit = bboptim(par, fn [, control] [, gr])
```

with corresponding `bboptim_box`, `bboptim_projected`, and `bboptim_linear`
variants.

## Nonlinear systems

```fortran
fit = sane(par, fn [, control])
fit = dfsane(par, fn [, control])
fit = bbsolve(par, fn [, control])
```

`sane_result` fields are `par`, `residual`, `fn_reduction`, `feval`, `iter`,
`convergence`, `method`, `m`, `nm_used`, and `message`.

SANE convergence codes match its R routine:

- 0: success
- 1: maximum iterations
- 2: function-evaluation error
- 3: maximum line-search reductions
- 4: anomalous iteration
- 5: lack of improvement

DF-SANE uses the R DF-SANE codes:

- 0: success
- 1: maximum iterations
- 2: stagnation
- 3: function-evaluation error
- 4: maximum line-search reductions
- 5: lack of improvement

## Projections

```fortran
call project_box(x, lower, upper, projected, ok)
call project_linear(x, a, b, meq, projected, ok)
```

For `project_linear`, constraints are stored by rows and interpreted as
`A*x >= b`; rows `1:meq` are equalities.

## Multistart

```fortran
ans = multistart_solve(starts, fn [, control])
ans = multistart_optimize(starts, fn [, control] [, gr])
ans = multistart_optimize_box(starts, fn, lower, upper [, control] [, gr])
ans = multistart_optimize_projected(starts, fn, project [, control] [, gr])
ans = multistart_optimize_linear(starts, fn, a, b, meq [, control] [, gr])
```

Each column `starts(:,k)` is one starting point.
