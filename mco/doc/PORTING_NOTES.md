# Porting notes

## Algorithm

The optimizer follows NSGA-II: bounded random initialization, constrained
nondominated sorting, crowding distance, binary tournament selection,
simulated-binary crossover, polynomial mutation, and elitist truncation of
the combined parent/offspring population.

The upstream C code stores a negative sum of violated constraints. The
Fortran code stores the equivalent positive violation magnitude. In both
interfaces, a constraint callback is feasible when each returned value is
nonnegative.

## Deliberate adaptations

- The intrinsic Fortran random generator is seeded deterministically rather
  than using R's RNG state.
- The hypervolume routine is a newly written exact recursive slicing method.
  It is suitable for modest fronts and dimensions, but it does not reproduce
  the performance of the upstream specialized dimension-sweep C algorithm.
- `generalized_spread` uses absolute deviations from the nearest-neighbor
  mean, the standard published form. The upstream R expression omits the
  absolute value.
- R vectorized callbacks and collections saved at several requested
  generations are not emulated. A result contains the final generation.
- Plotting and S3 dispatch are omitted.

## Numerical conventions

Every objective is minimized. Parameters and objective values are column-wise:
`par(nvar,npop)` and `value(nobj,npop)`. The default mutation probability is
`1/nvar`, matching conventional NSGA-II behavior when a negative option value
is retained.
