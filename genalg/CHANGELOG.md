# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `genalg` 0.2.1.
- Translate bounded floating-point `rbga()` genetic algorithm.
- Translate binary `rbga.bin()` genetic algorithm.
- Preserve rank-biased selection, one-point crossover, elitism, suggestions,
  generation histories, evaluation caching, and package-specific mutation.
- Preserve the historical binary stale-evaluation behavior by default, with a
  switch for corrected invalidation.
- Add standalone deterministic RNG and objective/monitor callback interfaces.
- Add five regression programs and two examples.
- Retain original R sources/manuals/package metadata under `original/`.
- Preserve GPL-2 licensing.
