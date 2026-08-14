# Algorithm notes

## Generalized PAVA

The PAVA implementation follows the upstream block-merging algorithm. Each
adjacent violation merges two blocks and recomputes the selected block
functional. Mean, median and fractile solvers are supplied directly. Repeated
measurements retain the upstream behavior of explicitly merging the response
and weight vectors in a block.

Primary ties are ordered by predictor and response/block value. Secondary ties
collapse equal predictors before PAVA. Tertiary ties preserve within-tie
response deviations and add the isotonic correction of the collapsed group.

## Active-set method

The upstream primal active-set iteration is retained:

1. identify currently active order constraints;
2. solve the equality-constrained loss problem;
3. remove an active constraint when a Lagrange multiplier violates dual
   feasibility;
4. otherwise line-search from the current feasible point to the
   equality-constrained proposal and activate the first blocking constraint.

Connected components of active pairwise constraints provide the reduced
parameterization. LS/L1/quantile/Chebyshev and Lp block minimizers are solved
directly. Smooth losses use a reduced-space BFGS implementation in place of
R's `optim(method="BFGS")`.

The nonsmooth L1 and quantile solvers intentionally retain the upstream
subgradient convention (`sign(0)=0`). As in the R implementation, KKT
stationarity diagnostics can therefore be conservative when the optimum lies
exactly on observed response values.

## Regression with restricted fitted values

`mregnn`, `mregnnM`, and `mregnnP` follow the upstream dual construction:
an orthonormal basis for the column space of X is formed, the dual problem is
solved by nonnegative least squares, and fitted values are reconstructed from
the projected solution. The supplied `nnls-fortran-v0.1.0` implementation is
vendored unchanged except for its location in this package.

## Portability

The code does not rely on Fortran short-circuit Boolean evaluation. Bounds
checking was used specifically to audit branches where the R source combines
index tests with logical conditions.
