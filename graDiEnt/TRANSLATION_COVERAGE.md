# Translation coverage

Upstream: `graDiEnt` 1.0.1.

## Directly translated

- `GetAlgoParams`
- `optim_SQGDE`
- `SQG_DE_bin_1_rand`
- `SQG_DE_bin_1_curr`
- `SQG_DE_bin_1_best`
- `Purify`

The three adaptation schemes preserve the upstream stochastic quasi-gradient,
`psi` self-scaling, binomial crossover, base-vector selection, uniform jitter,
and greedy acceptance rules.

The population update is synchronous: all proposals in one iteration are
computed from the same pre-update population and objective vector.  This
matches both the upstream sequential `lapply` path and its parallel
`parLapplyLB` path.

## R-only code omitted

- `doParallel` registration;
- `parallel::makeCluster`, `parLapplyLB`, and cluster shutdown;
- console `message()` progress output;
- R lists, matrices, S3/runtime plumbing;
- plotting in examples.

Parallel execution changes performance rather than the mathematical update, so
it is not a required dependency of the Fortran core.

## Source-fidelity details

- `rand` samples `2*n_diff+1` mutually exclusive population indices; the last
  is the DE base vector.
- `current` samples `2*n_diff` parents excluding the target particle.
- `best` samples `2*n_diff` parents excluding the current best particle and
  uses the best particle as the base vector.
- The R source draws `n_params` uniform jitter values even when crossover
  selects fewer coordinates.  Its subset replacement consumes the first
  selected-count draws.  The Fortran translation preserves this draw count and
  ordering.
- Failed/undefined stochastic gradient normalization leaves the target
  particle unchanged before objective re-evaluation, as in the source.
- Initialization rejects all non-finite objective values.  Proposal evaluation
  maps NaN to a worst value; negative infinity remains an improving value,
  matching R's `is.na` rather than `is.finite` check in the adaptation code.

## Defensive fixes to upstream edge cases

These changes avoid runtime failures while preserving the intended algorithm:

1. **Parent-count validation.**  Upstream checks only
   `n_diff <= n_particles/2`, but all strategies can require
   `2*n_diff+1 <= n_particles`.  The Fortran validator uses the stronger
   condition instead of allowing `sample(..., replace=FALSE)` to fail later.
2. **Trace allocation.**  Upstream allocates `floor(n_iter/thin)` trace slots,
   but its index logic can require `ceiling(n_iter/thin)` when the division is
   not exact.  Fortran allocates the required ceiling count.
3. **Convergence with thinning.**  The R convergence code indexes the thinned
   trace using `stop_check` as though it were an iteration count, which can
   underflow when `thin>1`.  Fortran maintains a rolling window of the actual
   last `stop_check` iteration objective vectors.
4. **Zero jitter.**  Upstream documentation says `jitter_size=0` disables
   jitter, while `GetAlgoParams` rejects zero.  Fortran accepts zero, following
   the documented behavior.
5. **Subset jitter assignment.**  R's full-length jitter RHS can emit a
   recycling/replacement warning when crossover selects a subset.  Fortran
   reproduces the effective first-draw assignment without emitting a warning.

## Numerical differences

R and Fortran use different random-number generators and normal variate
implementations, so identical integer seeds do not imply identical particle
trajectories.  With corresponding initial states/random draws, the translated
update equations are otherwise the same except for the defensive fixes above.
