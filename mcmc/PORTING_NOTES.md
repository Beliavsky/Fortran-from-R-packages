# Porting notes

## Scope

This is a native modern-Fortran translation of the computational core of
`mcmc` 0.9-8. The upstream package has no plotting API; its presentation layer
is principally R S3 objects, formula-like callback closures, and debug/list
assembly. Those R object mechanics are replaced by typed Fortran callbacks and
derived result types.

## Random-walk Metropolis

The `metrop` kernel preserves the upstream iteration order:

1. perform `nspac` random-walk Metropolis transitions;
2. evaluate the output function once;
3. accumulate `blen` such outputs into a batch mean;
4. repeat for `nbatch` batches.

Acceptance rates are reported overall and per batch. Optional debug arrays
retain current state, proposal, log Hastings ratio, generated normal vector,
uniform draw (or `-1` when not needed), and acceptance decision.

For full matrix scales the proposal is `state + scale %*% z`, matching the C
implementation's column-major traversal.

## Simulated/parallel tempering

The two upstream state representations are separated into two strongly typed
entry points:

- `temper_serial`: state is `[component, x...]`;
- `temper_parallel`: state is a `ncomp x dimension` matrix.

Each transition chooses with probability 1/2 between a within-component
random-walk update and a jump/swap update. Serial jumps include the neighbor
proposal asymmetry correction

`log(degree(i)) - log(degree(j))`.

Parallel swaps evaluate both crossed component densities exactly as upstream.
`acceptx` and `accepti` preserve component-specific acceptance accounting;
`-1` is used where R would report `NA` because an edge/update was never
attempted.

The neighbor matrix must be symmetric and every component must have at least
one neighbor.

## Morphometric transformations

`morph_transform` translates the isotropic radial maps used by `morph`:

- polynomial/exponential expansion controlled by `r` and `power`;
- subexponential expansion controlled by `b`;
- their composition;
- optional scalar or vector center.

The forward map, inverse map, and inverse-map log Jacobian are available as
type-bound procedures. General powers use Newton iteration for the monotone
radial inverse. `morph_metrop` evaluates the original-space density plus the
Jacobian in transformed coordinates and reports batches in original
coordinates.

## Initial sequence estimators

`initseq` is a direct translation of `src/initseq.c`:

- initial-positive paired autocovariances;
- monotone sequence via cumulative minima;
- convex sequence via the same pool-adjacent-violators algorithm;
- `var_pos`, `var_dec`, and `var_con` use `2*sum(Gamma)-gamma0`.

As upstream documents, these estimates can rarely be negative for chains that
are much too short.

## Overlapping batch means

`olbm` is a direct translation of `src/olbm.c`, including the normalization
for the estimated covariance of the sample mean and the `demean=.false.` mode
that treats the true mean as zero.

## RNG and reproducibility

The R package uses R's `unif_rand()` and `norm_rand()`. The Fortran port uses
Fortran's intrinsic uniform RNG plus Box-Muller normals. `set_mcmc_seed`
therefore gives reproducibility within this port but does **not** reproduce R's
`.Random.seed` bit-for-bit.

R's `initial.seed`/`final.seed` list elements are not reproduced because the
standard Fortran RNG has no portable serialized-state representation. A caller
that needs exact stream checkpointing should provide its own RNG layer in a
future extension.
