# Third-party licensing and provenance

## PWEV

The original R package declares GPL-3. The retained source is in `original/PWEV-master`.

## rumidas Fortran dependency

The attached translation is vendored in `vendor/rumidas-fortran` and declares GPL-3.0-only. It includes its own `LICENSE`, notices, upstream snapshot, and vendored `maxLik` dependency.

## rugarch Fortran dependency

The attached translation is vendored in `vendor/rugarch-modern-fortran` and declares GPL-3.0-only. It includes its own licensing and provenance files.

## WeightedEnsemble and metaheuristicOpt algorithms

The original PWEV workflow calls the GPL-3 R package `WeightedEnsemble`, which in turn uses the GPL-2-or-later `metaheuristicOpt` PSO routine. This port implements only the PSO path needed by PWEV and distributes the combined work under GPL-3.0-only.
