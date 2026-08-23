# Licensing

## rangen

The supplied R package declares:

```text
License: GPL-3
```

The translated Fortran project is therefore distributed under **GPL-3.0-only**. The full GPL v3 text is in `LICENSE`.

Upstream authors:

- Manos Papadakis — author, maintainer, copyright holder
- Michail Tsagris — contributor
- Omar Alzeley — contributor

## PCG32 lineage

The upstream `Random.h` states that its minimal PCG32 code is derived from M. E. O'Neill's PCG code and is licensed under Apache License 2.0. The Fortran implementation translates the recurrence and output permutation rather than copying C++ syntax. Apache-2.0 code is compatible with distribution of the combined work under GPL-3.0.

## zigg

The upstream package depends on the external R package `zigg` for its Ziggurat normal generator. No `zigg` source is bundled in the supplied archive and no `zigg` source is redistributed here. The Fortran port uses its own Box-Muller normal implementation instead.
