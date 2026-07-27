# Porting notes

## Representation

The R package appends probabilities as the last column of `dat`. The Fortran
interfaces accept `returns(:, :)` and `probabilities(:)` separately. This avoids
mixing observations and metadata in one array.

R lists are represented by typed result structures. Fortran functions return
full precision and do not round numerical fields to a tolerance-dependent number
of decimal places.

## LP solver replacement

The original Benders master problem calls `Rsymphony_solve_LP`. The Fortran port
contains a two-phase tableau simplex implementation for

```text
minimize c^T x
subject to A x <= b
           x >= 0.
```

It supports negative right-hand sides through surplus and artificial variables,
uses a Phase-I feasibility objective, removes artificial variables, and reports
infeasible or unbounded problems explicitly.

## Benders decomposition

The scenario cut generation follows the original centered-return algorithm. An
important ordering detail is preserved: the lower bound in each Benders
iteration is computed using the active scenario set that generated the current
master cut, while the upper bound uses the newly classified active set.

## Projection

The original package applies the Zhao-Li path-following equations. A direct
translation is available as `zi_projection` and the residual equations are
available as `f_func`.

For the high-level routine, the default is a more robust two-stage method:

1. solve the complete risk LP with the native simplex solver;
2. minimize the weighted squared distance to the benchmark subject to the LP
   constraints and the optimal-risk bound using dense ADMM.

This method correctly recovers a benchmark that already lies on a degenerate
optimal face. `use_zi=.true.` selects the translated path-following method.

## Constraint-sign correction

The public R documentation states:

```text
Aconstr * theta <= bconstr
```

`PortfolioOptimProjection`, however, passes `Aconstr` into an internal system
whose limiting equation represents `A x >= b`. The Fortran high-level routine
follows the documented less-than-or-equal convention by default. Set
`upstream_constraint_sign=.true.` to reproduce the internal R sign convention.

## Lower-bound benchmark correction

The projection LP solves for shifted weights `x = theta - LB`. The upstream code
passes the unshifted benchmark directly as `xhat`. The Fortran default projects
toward `benchmark - LB`, which is consistent with the optimization variables.
Set `upstream_benchmark_shift=.true.` for the original behavior.

## Tail-risk boundary correction

The upstream `.RISK_post` expression forms the range `(index + 1):n`. In R, this
can create a descending or out-of-range sequence when the VaR index is the final
scenario. The Fortran implementation handles an empty strict tail explicitly.

## Target-return adjustment

When an upper-bound-only calculation shows that the requested return is too
high, the upstream code subtracts one twentieth of the original target from the
bound-based maximum. Large targets can therefore produce a severely negative
replacement target. The Fortran port instead moves the target just below the
computed bound-based maximum by a tolerance-scaled amount.

## Probability handling

Nonnegative probabilities are normalized internally. Zero total probability and
negative probabilities are rejected.

## Numerical limitations

The simplex implementation is dense. The Benders master remains small even for
large scenario counts, but the projection formulation contains one auxiliary
variable and constraint per scenario and is intended for moderate sample sizes.
The ADMM second stage also uses dense normal equations.
