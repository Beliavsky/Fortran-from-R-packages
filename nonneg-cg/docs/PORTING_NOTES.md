# Porting notes

## Directly retained algorithm

The search direction, modified Polak-Ribiere-Polyak coefficients, projection at
the nonnegative boundary, maximum feasible step, and squared-step sufficient
decrease test follow `src/nonnegcg.c` from the supplied package.

## Intentional Fortran API additions

- Explicit callback interfaces and optional polymorphic user data.
- A monitor callback that can cancel a run.
- Validation of dimensions, controls, feasibility, and finite callback output.
- `objective_calls`, which gives the true objective-call count.

## Original evaluation-count behavior

The original C code initializes `nfeval` to one and increments it only after a
line-search trial is rejected. Accepted line-search evaluations are therefore
not included. This port preserves that value in `result%nfeval` for regression
compatibility and separately records every objective call in
`result%objective_calls`.

For the documented Rosenbrock problem, the original C kernel reports 372 while
actually calling the objective 419 times.

## Parallelism

The original source contains optional OpenMP loops and attempts to change BLAS
thread counts. The numerical method is inexpensive per vector element and the
Fortran port is dependency-free, so thread-runtime manipulation is omitted.
Applications remain free to parallelize expensive objective and gradient
callbacks themselves.
