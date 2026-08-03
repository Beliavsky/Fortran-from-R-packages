# Porting notes

## Numerical equivalence

The upstream `matrixFDM` function stores a fixed set of sparse diagonals. This
translation constructs the same differential operator row by row from
one-dimensional nonuniform derivative stencils:

- central first and second differences at interior nodes,
- inward first-order differences at boundaries,
- no boundary second-difference term,
- tensor products of first derivatives for correlation terms,
- subtraction of the risk-free discount rate on the diagonal.

This formulation is easier to audit in Fortran and supports different node
counts for each asset without reproducing R `expand.grid` objects.

## Sparse solver

R delegates linear systems to the `Matrix` package. The Fortran port uses a
native CSR matrix and Jacobi-preconditioned BiCGSTAB. Its convergence test uses
a diagonally scaled residual. Scaling is important for American options because
the penalty diagonal can be many orders of magnitude larger than continuation-
region coefficients.

## American exercise

The penalty equation and stopping conditions follow the upstream algorithm.
Each step stops when the active set is unchanged, the relative iterate change
is below tolerance, or the configured penalty-iteration limit is reached.

## Deliberate interface changes

- R lists are replaced by derived types.
- Option value arrays are always stored as rank-one classical state vectors.
- Error conditions are returned through `status_type` rather than R exceptions.
- Progress messages are not emitted from library routines.
- Plotting is omitted.

## Computational limits

The method is not an ADI scheme. The total node count is the product of all
per-asset node counts, and mixed derivatives increase the sparse bandwidth.
The port is intended mainly for low-dimensional multi-asset problems, as is the
upstream package.
