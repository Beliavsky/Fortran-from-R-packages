# Algorithm notes

## Upstream algorithm retained

For each offspring, two tournaments are sampled without replacement from the
current population. Within a tournament, candidates receive the same average
rank weights as `rank(-fitness)` in R, so lower objective values receive larger
selection weights. A parent is sampled with probability proportional to that
rank.

Crossover forms the ordered union of the two parent subsets and uniformly
samples `k` distinct elements from it.

Each of the `k` child positions mutates independently with probability
`mutprob`. If `r` positions mutate, `r` distinct indices are sampled from
`1:n` excluding every index currently in the child. This guarantees fixed
cardinality and no duplicate indices.

When elitism is enabled, the `keepbest` best members of the old population and
the `popsize-keepbest` best offspring form the next generation. Random tie
keys reproduce the upstream `ties.method="random"` intent.

## `mutfrac`

The upstream conversion is preserved exactly:

`mutprob = 1 - (1 - mutfrac)^(1/k)`.

Thus `mutfrac` controls the expected fraction of offspring experiencing at
least one mutation.

## RNG

The R package delegates randomness to R's global RNG. This port uses the
minimal-standard Park-Miller generator with Schrage arithmetic, avoiding
integer-overflow assumptions and making the Fortran seed reproducible across
compilers. Exact seeded R/Fortran trajectories are therefore not expected to
match, but the probability laws and GA rules do.

## Parallel evaluation

Upstream parallelization affects only objective evaluation and is implemented
with R's `parallel` and `bigmemory` packages. It is not part of the genetic
operator mathematics and is intentionally outside this translation.
