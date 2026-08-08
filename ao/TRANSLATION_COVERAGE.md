# Translation coverage

Upstream: `ao` 1.2.3.

## Translated package-owned computation

- `ao()` single-process iteration logic.
- Sequential, random, none, and custom partitions.
- `generate_random_partition()` including shuffling, Bernoulli block starts,
  and forcing the requested minimum block count.
- `split_by_target()` numeric splitting semantics (`split_estimate`).
- `Process` numerical state/history and best-point selection.
- Value, parameter, iteration, time, and subproblem-error stopping logic.
- `merge_results()` numerical best-process selection, represented by
  `ao_multi_result`.
- Minimization/maximization, bounds, gradients, Hessians, and custom parameter
  norm callbacks.

## R/runtime infrastructure intentionally omitted

- R6 classes and active bindings;
- formula/function-argument introspection supplied by `optimizeR::Objective`;
- named target argument reconstruction and arbitrary `...` forwarding;
- `checkmate`, `cli`, `oeli`, and `progressr` validation/UI;
- `future.apply` parallel scheduling;
- package startup messages and documentation/graphics.

The numeric `npar` splitting helper is retained, but Fortran does not emulate
R's named argument lists.

## External base optimizers

`ao` does not implement its base optimization algorithms. It delegates to
`optimizeR`, with `stats::optim(method="L-BFGS-B", maxit=10)` by default.
Therefore an exact translation of `ao` cannot truthfully label a newly written
L-BFGS-B as upstream `ao` code.

This standalone release supplies bounded BFGS, Nelder-Mead, and Newton block
solvers so the AO algorithm is directly usable. They are compatibility
implementations. The AO orchestration, block construction, stopping logic, and
history semantics are the translated package-owned computation.

## Timing

R records the `seconds` reported by each external optimizer. The Fortran port
records `cpu_time` around each block solve. Time-limit checks remain between
subproblems, as in upstream.

## Random-number streams

The partition algorithm matches the upstream sampling logic conceptually, but
Fortran's intrinsic RNG is not R's Mersenne-Twister stream. Equal integer seeds
do not imply identical R/Fortran partitions.
