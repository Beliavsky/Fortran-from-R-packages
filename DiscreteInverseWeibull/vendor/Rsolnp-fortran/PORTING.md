# Porting notes

## Scope

The nine exported R entry points are represented in the Fortran API:
`solnp`, `csolnp`, `csolnp_ms`, `gosolnp`, `startpars`, `kkt_diagnose`,
`solnp_standardize_problem`, `solnp_problem_suite`, and
`solnp_problems_table`.

R printing, S3 classes, formula/list dispatch, `future.apply` parallelism,
truncated-normal random starts, vignettes, and documentation-generation code
were omitted. Original R and C++ computational sources remain under `original/`.

## Solver adaptation

The upstream package has legacy R and newer C++ paths. The Fortran translation
uses one self-contained augmented-Lagrangian implementation for both `solnp`
and `csolnp`. Major iterations update constraint multipliers and penalties;
each subproblem is minimized by a projected BFGS method with Armijo
backtracking. Two-sided nonlinear inequalities are represented internally by
bounded slack variables.

The translated solver estimates final multipliers by a regularized
least-squares KKT system using equality and active-inequality Jacobians. This
prevents stale penalty multipliers from making stationarity diagnostics appear
poor after a correct primal solution has been found.

## Derivatives

Associated analytic gradient/Jacobian callbacks are preferred. Missing
callbacks use central finite differences with scale `control%delta`. This is a
compact derivative layer adapted from the design used in the earlier
`numDeriv-fortran` translation; it is not the complete Richardson/complex-step
API of that project.

## Multistart

The R package can use random sampling and parallel workers. This translation
uses deterministic Halton points, optional feasibility improvement, and
sequential solves. `seed` offsets the low-discrepancy sequence rather than
seeding a global random-number generator. Fixed parameters are represented by
equal lower and upper bounds.

## Standard form

The R helper returns nested lists describing transformed functions and bounds.
Fortran returns a `solnp_problem` with explicit standard-form metadata. Equality
targets are shifted to zero, and each finite lower/upper nonlinear inequality
bound contributes a one-sided constraint entry.

## Benchmark suite

The complete 77-row registry is available, and all original benchmark source is
retained. Eighteen representative benchmark definitions have executable
Fortran callbacks. Requesting an unimplemented registry entry returns
`solnp_invalid_problem` instead of silently substituting another problem.

The GARCH benchmark uses a portable Park-Miller/Box-Muller generated series, so
its exact sample and objective differ from R's RNG stream. The mathematical
likelihood and parameter constraints are retained.

## Result differences

The result type is explicit rather than an R list. The Hessian component is the
final BFGS approximation to the augmented objective, not an exact Hessian of the
original objective. Timing uses `cpu_time`; values can differ from wall-clock
measurements in R.
