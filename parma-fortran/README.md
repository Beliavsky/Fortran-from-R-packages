# parma-fortran

A modern Fortran 2018 and FPM translation of the computational parts of the R
package `parma` 1.7, Portfolio Allocation and Risk Management Applications.
The supplied original package is retained unmodified under
`original/parma-1.7`.

## Build

```text
fpm build
fpm test
fpm run parma_demo
fpm run --example risk_measures
fpm run --example efficient_frontier
```

The library has no external runtime dependencies. It uses allocatable arrays,
derived types, procedure callbacks, explicit interfaces, and `implicit none`.

## Main API

```fortran
use parma

type(parma_spec) :: spec
type(parma_port) :: solution
type(parma_options) :: options
real(dp) :: returns(1000, 5), lb(5), ub(5)
integer :: info

lb = 0.0_dp
ub = 1.0_dp
call parmaspec(spec, data=returns, risk=risk_cvar, &
   objective=solve_min_risk, lb=lb, ub=ub, alpha=0.05_dp, info=info)
call parmasolve(spec, solution, options)
print *, solution%weights, solution%risk, solution%reward
```

`parmafrontier` solves a sequence of target-return problems. `parmautility2`
and `parmautility4` implement the second- and fourth-moment CARA
approximations. The umbrella module also exports:

- MAD, variance, minimax, CVaR, CDaR, LPM, UPM, and Rachev measures.
- Scenario transformations and benchmark-relative covariance risk.
- Budget, leverage, target, turnover, buy/sell turnover, variance, linear, and
  cardinality constraints.
- Full-covariance CMA-ES with bound reflection.
- Convex box/budget QP, canonical continuous simplex LP, exact small binary
  optimization, and a general SOCP logarithmic-barrier solver.
- Lag and business-day routines, matrix helpers, entropy, simulated feasible
  weights, and lightweight result extractors.

## Risk-measure compatibility

The original public R `fun.lpm` expression uses the opposite tail from the NLP
optimizer. The Fortran API defaults to the conventional downside definition,
`max(threshold-return, 0)`. Set `spec%lpm_legacy = .true.` or pass
`legacy=.true.` to reproduce the public R expression. A threshold of `999`
centers portfolio returns before evaluating the partial moment, as documented
by the original package.

## Scope

This is a computational translation, not an R compatibility layer. R S4
classes, plotting, `xts` integration, parallel-cluster orchestration, and calls
to R packages such as `nloptr`, `Rglpk`, and `quadprog` are not dependencies.
Native Fortran solvers and typed data structures replace them. See
`COVERAGE.md` and `PORTING_NOTES.md` for exact mappings and limitations.

## License

GPL-3.0-or-later. See `LICENSE`, `COPYRIGHTS`, `NOTICE`, and the SPDX headers in all
Fortran sources.
