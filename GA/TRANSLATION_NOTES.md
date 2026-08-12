# Translation notes

## Source version

This release was translated from the attached `GA` 3.2.5 source tree.

## Architectural mapping

### `R/ga.R`

The main generation loop maps to `ga_core.f90`:

- population initialization and suggestions;
- lazy fitness evaluation after crossover/mutation;
- six-column GA fitness summary;
- `garun` and `maxFitness` stopping;
- selection, mating, crossover, mutation, and elitism;
- final best solution and fitness;
- optional best-solution history.

Fortran uses explicit result/control derived types instead of an R S4 object.

### `src/genope.cpp` / `R/genope.R`

The package's native genetic operators map to `ga_operators.f90`. All exposed
selection families, binary operators, real-valued operators, permutation
operators, adaptive mutation probability, and the DE update are represented.

A few low-level Rcpp helper functions such as `which_asR`, vector concatenation,
and set intersection/difference are implementation details and are expressed
idiomatically with Fortran loops instead of being public APIs.

### `R/gade.R`

`de_real` implements the package DE/rand/1/bin update:

1. sample three distinct population members;
2. form `v = x_r1 + F*(x_r2-x_r3)`;
3. binomial crossover with one forced coordinate;
4. randomly repair coordinates falling outside bounds;
5. greedily retain the candidate only when fitness improves.

The original implementation allows the target index itself to appear among the
three sampled indices, and this translation retains that behavior.

### `R/gaIslands.R`

`ga_islands_*` implements serial coarse-grained islands. Migration sources are
kept immutable during a migration epoch, matching the R code's use of the
completed `GAs` list while writing to a separate next-population list.

## Numerical/semantic details deliberately preserved

- GA maximizes fitness, just like the R package. Minimize an objective by
  returning its negative.
- The default real selection is linear fitness scaling.
- The default real crossover is local arithmetic crossover.
- The default real mutation is uniform random single-coordinate mutation.
- The default permutation operators are order crossover and inversion mutation.
- `garun` follows the package literally: it counts all recorded best-fitness
  values within `sqrt(epsilon)` of the historical maximum.
- `ga_pmutation` uses exponential decay with defaults `p0=0.5`, `p=0.01`.
- Laplace crossover defaults to `a=0`, `b=0.15`.
- BLX defaults to alpha `0.5`.
- Power mutation defaults to exponent `10`.
- PBX chooses positions as the unique set from `n` draws with replacement.
- PMX completion preserves donor-order set difference, as in the Rcpp code.
- Displacement mutation reinserts after at least one remaining element, matching
  the original implementation.

## Consolidations and differences

### Random-number generator

R's RNG API is replaced by the Fortran intrinsic RNG with deterministic seeding.
The distributions are equivalent but the exact streams differ.

### Elitism

The original R code avoids inserting duplicate elite rows when enough unique
rows exist. v0.1.0 retains the best `elitism` rows directly. This can only alter
behavior when the elite prefix contains duplicate chromosomes; best-fitness
preservation is unchanged.

### Fitness reuse after generic single-point crossover

The Rcpp generic crossover can reuse parent fitness when the crossover point is
at either extreme and therefore a child is unchanged. The Fortran generation
engine conservatively marks crossed individuals for reevaluation. The resulting
population/fitness values are unchanged, but evaluation counts can be slightly
higher.

### Multiple tied solutions

R's final S4 object may retain several distinct chromosomes tied for the best
fitness. `ga_real_result%solution` and `ga_int_result%solution` currently retain
the first best chromosome; complete final populations and fitness vectors are
available to recover all ties.

### Hybrid local optimization

`ga(..., optim=TRUE)` delegates to R's external `stats::optim`. Since `optim` is
not an algorithm implemented by GA itself, it is not duplicated here. Native
Fortran users can run a local optimizer on `result%solution` or embed one in a
higher-level workflow.

### Parallel execution

R `foreach`/cluster management is omitted. The island implementation is serial
and can be parallelized by applications at the island level if desired.

## Validation

Tests cover:

- binary/decimal and binary/Gray conversions;
- repair and reflection helpers;
- adaptive mutation probability;
- all six selection strategies;
- every translated crossover and mutation family, including permutation
  validity invariants;
- real-valued GA optimization and supplied-suggestion/max-fitness stopping;
- binary OneMax optimization;
- permutation optimization;
- Differential Evolution;
- real, binary, and permutation island GAs.

The package is tested both with runtime bounds checking and with `-O2`.
