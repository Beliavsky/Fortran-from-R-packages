# ao-fortran

Modern Fortran/FPM translation of the computational core of the R package
`ao` 1.2.3 (Alternating Optimization).

## Scope

The translation implements the package-owned alternating-optimization logic:

- sequential coordinate/block optimization;
- a newly generated random partition on every AO iteration;
- no partition (joint optimization);
- arbitrary custom blocks, including overlapping blocks;
- minimization and maximization;
- scalar or vector parameter bounds;
- optional analytic gradients and Hessians;
- value, parameter, iteration, and elapsed-CPU stopping criteria;
- the same `tolerance_history` comparison convention as the R `Process` class;
- detailed per-block history (iteration, value, full parameter vector, active
  block flags, and subproblem CPU time);
- returning the best point seen in the process, not merely the final point;
- multiple-process Cartesian products over starting values, partition modes,
  and base-optimizer choices;
- random-partition generation and numeric target splitting helpers.

The public code works directly with numeric vectors rather than R function
argument names, model objects, R6 classes, data frames, or `...` arguments.

## Base optimizer boundary

The R package delegates every block subproblem to the separate `optimizeR`
package. Its default is `stats::optim(method="L-BFGS-B", maxit=10)`. Those
optimizers are not code owned by `ao`.

To make this FPM package standalone, three compatible block solvers are
provided:

- `AO_BASE_BFGS` (default): bounded projected quasi-Newton with Armijo search;
- `AO_BASE_NELDER_MEAD`: bounded derivative-free simplex search;
- `AO_BASE_NEWTON`: analytic-Hessian Newton with ridge rescue and BFGS
  fallback if no Hessian is supplied.

This standalone base-solver layer is not claimed to be a line-for-line port of
R's `stats::optim` or `stats::nlm`. The AO orchestration surrounding it is the
translated package algorithm. The interface is deliberately separated so a
translated external optimizer can be substituted later.

## Example

```fortran
use ao
implicit none
type(ao_result) :: result
real(dp) :: initial(2)
initial = [2.0_dp, 2.0_dp]
call ao_optimize(rosenbrock, initial, result)
```

See `example/` and `API.md` for complete programs and callback signatures.

## Build

```text
fpm build
fpm test
fpm run --example rosenbrock
```

A strict GNU Fortran test script is available as `scripts/test_gfortran.sh`
and `scripts/test_gfortran.bat`.

## License

The original `ao` package is GPL-3.0-only. The full upstream source supplied by
the user is preserved under `original/ao-main/`, and the upstream `LICENSE` is
retained at the package root.
