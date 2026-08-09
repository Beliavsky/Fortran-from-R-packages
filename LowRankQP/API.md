# LowRankQP-fortran API

## Modules

### `lowrankqp_kinds`

Exports `dp = kind(1.0d0)`.

### `lowrankqp`

The main public routine is:

```fortran
call solve_low_rank_qp(v, d, a, b, u, result, options)
```

It solves

```text
min  d^T alpha + 1/2 alpha^T H alpha
s.t. A alpha = b
     0 <= alpha <= u
```

where `v` has shape `(n,m)` and

- if `m == n`, `H = v`, matching the upstream package;
- otherwise, `H = v * transpose(v)`.

Unlike the R wrapper, the Fortran constraint matrix is passed naturally as
`a(p,n)`, so each row of `a` is one equality constraint. The R wrapper passes
`t(Amat)` to the C core internally; this Fortran interface removes that
transposition artifact.

### Methods

```fortran
LRQP_LU
LRQP_CHOL
LRQP_SMW
LRQP_PFCF
```

or use `lowrankqp_method("CHOL")`, etc.

`LU` and `CHOL` factor the full `H + D` system. `SMW` uses the
Sherman-Morrison-Woodbury identity. `PFCF` is the product-form Cholesky
factorization translated from the package source.

### Options

```fortran
type(lowrankqp_options) :: opt
opt%method   = LRQP_PFCF
opt%max_iter = 200
opt%tol      = 1.0e-8_dp
opt%verbose  = .false.
```

The defaults match the R `LowRankQP()` defaults.

### Result

```fortran
type(lowrankqp_result) :: res
```

Fields:

- `alpha` -- primal solution.
- `beta` -- equality-constraint dual variables.
- `xi` -- upper-bound dual variables.
- `zeta` -- lower-bound dual variables.
- `iterations` -- number of interior-point statistics/termination passes.
- `status` -- `0` converged; `1` iteration limit; negative values are invalid input; positive values above 1 are factorization/solve failures.
- `converged` -- logical convergence flag.
- `primal_feasibility`.
- `dual_feasibility`.
- `complementarity`.
- `duality_gap`.
- `termination` -- the package's equation-(12)-style stopping quantity.

### Objective utility

```fortran
f = lowrankqp_objective(v, d, alpha)
```

uses the same square-vs-low-rank Hessian convention as the optimizer.
