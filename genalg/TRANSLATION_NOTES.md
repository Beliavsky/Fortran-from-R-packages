# Translation notes

## Scope

`genalg` 0.2.1 contains no compiled C, C++, or Fortran source. Its algorithms
are implemented in two R files:

- `R/rbga.R`
- `R/rbga.bin.R`

Those computational routines were translated. The following R infrastructure
was intentionally omitted:

- `plot.rbga` graphics;
- `summary.rbga` presentation logic;
- S3 class creation and dispatch;
- console-only `showSettings` formatting;
- R RNG state and R's `sample`, `runif`, and `dnorm` implementations.

The optional R monitor function is represented as an explicit Fortran callback.

## Representation

Fortran populations are stored as

```text
population(pop_size, chromosome_size)
```

so one chromosome occupies a row, matching the conceptual R representation.

The real and binary optimizers have separate result/control types because this
keeps callback interfaces type-safe.

## Parent selection

The R implementation sorts chromosomes from lowest to highest evaluation and
assigns parent weights

```text
dnorm(1:popSize, mean=0, sd=popSize/3)
```

The normalizing constant is common to all ranks, so the Fortran implementation
uses the mathematically equivalent weights

```text
exp(-0.5*(rank/(pop_size/3))**2)
```

and samples two distinct parents sequentially without replacement.

## Crossover and cached evaluations

The crossover point is uniform on `0:nvar`. Endpoint values copy a complete
parent and therefore copy its already-known fitness. Interior points produce a
new chromosome and require a new objective evaluation. This caching rule is
preserved.

For `nvar=1`, the original code does not crossover. Instead, it samples
chromosomes without replacement from the *sorted* population. The translation
preserves that branch explicitly.

## Real mutation

The exact original expression is retained:

```text
mutationVal = stringMax[var] - stringMin[var]*0.67
```

The mutation magnitude is multiplied by `(iters-iter)/iters` and a random sign.
An out-of-domain result is replaced by a fresh uniform draw from the gene's
bounds. Any real mutation invalidates the chromosome's cached evaluation.

## Binary initialization and mutation

The R code samples from `c(rep(0, zeroToOneRatio), 1)`. The Fortran port treats
`int(zero_to_one_ratio)` as the number of zero entries in this sampling pool,
which matches the intended count semantics of `rep(..., times=...)` for normal
nonnegative inputs.

Initial random binary chromosomes are repeatedly generated until at least one
bit is set. The same check is made for an interior-crossover child. As in the R
code, later mutation itself does not enforce the nonzero condition.

## Historical binary fitness-cache quirk

The original `rbga.bin()` mutates genes but never executes the real-valued
routine's equivalent of

```text
evalVals[object] = NA
```

Thus an endpoint-crossover child may be mutated while retaining its parent's
cached objective value. This behavior is part of the executable package code,
so v0.1.0 preserves it by default through
`legacy_binary_eval_cache=.true.`. Setting the flag false gives the logically
corrected behavior and invalidates any mutated chromosome.

## Random numbers

Exact R-seed reproducibility is not attempted. A deterministic standalone
Park-Miller generator provides uniform draws. This makes the library portable
and reproducible without linking against R while preserving the intended
sampling distributions.

## Enhancements that do not alter the algorithm

The result types expose `best_chromosome`, `best_value`, `nfe`, and
`mutation_count` as convenience diagnostics. The original R object already
contains enough information to derive the best chromosome/value, but not the
explicit objective-call count.
