# Computational coverage

## Public R entry points

| Original entry point | Fortran equivalent | Status |
|---|---|---|
| `parmaspec` | `parmaspec`, `type(parma_spec)` | Implemented |
| `parmasolve` | `parmasolve`, `type(parma_port)` | Implemented with native solvers |
| `parmafrontier` | `parmafrontier` | Implemented |
| `parmautility`, `parma2CARA`, `parma4CARA` | `parmautility2`, `parmautility4`, CARA value/gradient routines | Implemented |
| `riskfun` | `riskfun`, `risk_value` | Implemented |
| `cmaes` | `cmaes_minimize` | Full-covariance implementation |
| `Socp` | `socp_solve` | Native primal logarithmic-barrier implementation |
| `SocpControl` | arguments to `socp_solve` | Implemented in typed form |
| `weights`, `parmarisk`, `parmareward`, `parmastatus` | `parmaweights`, `parmarisk`, `parmareward`, `parmastatus` | Implemented |
| `checkarbitrage` | `checkarbitrage` | Implemented |
| exported turnover/variance constraints and Jacobians | typed residual and penalty routines in `parma_constraints` | Computational equivalents implemented |

## Risk and utility layer

Implemented:

- Shannon and cross entropy.
- Mean absolute deviation and population expected variance.
- Quadratic covariance risk and benchmark-relative covariance risk.
- Minimax loss, CVaR, CDaR, lower and upper partial moments.
- Rachev ratio and scenario-level MAD, variance, CVaR, and LPM transforms.
- Cumulative returns, drawdown/run-up series, empirical quantiles.
- Smooth maximum and absolute-value approximations and derivatives.
- Second- and fourth-moment CARA objectives and analytic gradients.

## Constraint layer

Implemented:

- Box bounds, budget, leverage, target equality or inequality.
- General linear equalities and two-sided inequalities.
- Total turnover and separate buy/sell turnover.
- Variance caps and maximum-position cardinality.
- Constraint feasibility checks, squared penalties, and weight repair.

The R package accepts arbitrary R functions and Jacobians. Fortran callers can
express arbitrary objectives through the `objective_callback` interface and
`cmaes_minimize`; the high-level `parma_spec` deliberately keeps constraints
typed rather than storing heterogeneous R closures.

## Solver layer

- `parma_cmaes`: full covariance adaptation, cumulative step-size adaptation,
  bound reflection, deterministic seeding, and stopping tests.
- `parma_qp`: convex quadratic programming with box, budget, and optional target
  constraints using projected gradient and line search.
- `parma_lp`: simplex for canonical continuous problems
  `max/min c'x`, `A x <= b`, `b >= 0`, `x >= 0`.
- `parma_milp`: exact enumeration for binary models up to a configurable small
  dimension; high-level portfolio cardinality uses a repair/penalty method.
- `parma_socp`: general constraints
  `||A_i x+b_i|| <= c_i'x+d_i` through phase-one feasibility and a Newton
  logarithmic barrier.
- `parmasolve`: deterministic projected variance optimization when applicable;
  otherwise constrained CMA-ES plus coordinate polishing.

## Time-series and numerical helpers

Implemented business-day sequences, weekday calculations, month-end indices,
vector/matrix lagging, covariance, symmetric eigensystems and square roots,
positive-definite adjustment, linear solves, projections, matrix triangles,
condition numbers, percentiles, bound reflection, and feasible random weights.

## Deliberately not translated as runtime functionality

- R S4 method dispatch, slots, formula parsing, `xts`/`timeSeries` classes.
- Plotting and PostScript output in the upstream test script.
- Parallel R cluster orchestration.
- The serialized `etfdata.rda` object. It remains in the retained original tree.
- Exact argument-for-argument wrappers around external R packages.

The original LP/MILP/QP/SOCP setup routines construct large solver-specific
slack-variable matrices. Their portfolio objectives and constraints are
available through the typed high-level API, but those R-specific intermediate
matrix layouts are not reproduced verbatim.
