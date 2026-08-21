# Translation notes

## Source baseline

- R package: DEoptim 2.2-8
- Package license: GPL (>= 2)
- Main numerical source: `src/de4_0.c`
- Population evaluator: `src/evaluate.c`
- R control/mapping logic: `R/DEoptim.R`

Reference copies are retained under `original/`.

## C/R to Fortran mapping

The C engine stores an `NP x D` R matrix as a one-dimensional column-major
array. The Fortran implementation uses an explicit `(NP,D)` allocatable array,
which preserves the same population/member interpretation without pointer
arithmetic.

R callable objects are replaced by explicit Fortran procedure interfaces:

- `de_objective(x)` returns one scalar objective value.
- optional `de_map(x)` modifies one population member in place.

The R wrapper applies `fnMap` row by row, detects duplicate mapped rows, and
retries duplicate rows with random candidates up to five times. That behavior
is implemented in `map_population`.

## Differential-evolution details retained

The mutation equations for strategies 1 through 6 are translations of the
active switch in `de4_0.c`. Crossover starts at a random coordinate, mutates at
least one coordinate, and proceeds cyclically while a uniform draw is below
`CR`. Out-of-range trial coordinates are redrawn uniformly within the relevant
bound interval, matching the C source.

Strategy 6 sorts the previous generation by objective value and selects a
p-best member from the best `round(p*NP)` members, with a minimum of two.

When `c > 0`, the C source draws `CR` from a normal distribution and `F` from a
Cauchy distribution. The translation deliberately retains two source-level
quirks for compatibility:

1. the successful-parameter accumulators are not reset each generation;
2. because the C engine stores only scalar `d_cross` and `d_weight`, the final
   member's adapted values are credited to all successful trials during that
   generation.

These details differ from some textbook JADE implementations but reflect the
actual DEoptim 2.2-8 computational source.

## Stopping and histories

The loop condition and tolerance test follow the C engine. `bestmemit` and
`bestvalit` record the best point at the beginning of each generation, before
that generation's trial population is selected. The returned `bestmem` and
`bestval` include improvements from the final completed generation.

`nfeval` counts `NP` objective calls for the initial population and `NP` for
each completed generation. Mapping calls are not counted as objective calls.

## Differences from R

- R's exact random-number stream is not reproduced. This library uses a
  standalone deterministic RNG so it has no R runtime dependency.
- R warnings caused by control normalization are silently normalized in the
  Fortran `de_control`; invalid structural inputs return a status/message.
- `trace` defaults to zero in the Fortran type to keep library use quiet; set a
  positive integer for DEoptim-style periodic progress output.
- R-specific parallel dispatch is outside the numerical engine and omitted.
- The result uses allocatable arrays rather than an S3 list.
- Population-history arrays are trimmed to the generations/populations actually
  returned rather than preserving unused R allocation slots.

## `bs`

The DEoptim help page describes best-of-parent-and-child selection, but the
active C code executes `error("bs = TRUE not currently supported")`. The
Fortran translation therefore reports `de_unsupported` if `control%bs` is true.
