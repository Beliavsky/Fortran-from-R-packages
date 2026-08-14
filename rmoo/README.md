# rmoo-fortran v0.1.0

Modern Fortran/FPM translation of the computational code in the R package
`rmoo` 0.3.2.

The library implements multi- and many-objective evolutionary optimization for
real, binary, integer/discrete, and permutation decision variables.  All source
is free-form Fortran 2018 and the public convenience module is `rmoo`.

## Implemented algorithms

- NSGA-I, including Pareto ranking and fitness sharing
- NSGA-II fast nondominated sorting, crowding distance, tournament selection,
  and elitist environmental selection
- NSGA-III Das-Dennis reference directions, ideal/worst/extreme/nadir updates,
  reference-line association, niche counts, and niching
- R-NSGA-II normalized preference distance and epsilon-based preference
  truncation
- generational distance, inverted GD, GD+, and IGD+
- reference-point generation, scaling, and multi-layer merging
- real SBX crossover and polynomial/random mutation
- single-point, uniform, and HUX crossover
- permutation ordered crossover and inversion mutation
- binary and discrete integer mutation
- optional initial suggestions and deterministic seeding

The high-level routines are:

```fortran
call rmoo_optimize_real(...)
call rmoo_optimize_binary(...)
call rmoo_optimize_integer(...)
call rmoo_optimize_permutation(...)
```

Select the algorithm with `ALG_NSGA1`, `ALG_NSGA2`, `ALG_NSGA3`, or
`ALG_RNSGA2`.  Objective vectors are minimized, matching rmoo's
nondominated-sorting convention.

## GA dependency

The supplied `GA-fortran-v0.1.0` translation was used.  Only the modules needed
by rmoo (kind/RNG/utilities and GA genetic operators) are vendored under
`src/vendor_ga`.  The scalar GA/DE/island optimizers are not duplicated because
rmoo's multiobjective evolutionary loop is implemented natively here.

## Build

```text
fpm build
fpm test
fpm run --example zdt1_example
```

The manifest explicitly uses free source form and disables implicit typing and
implicit external procedures.

## Intentionally omitted R infrastructure

Plotting, S4 classes, print/summary/progress methods, `foreach` parallel
clusters, Plotly/ggplot2 integration, and R runtime validation/dispatch are not
translated.  Fortran derived result types expose population, objective values,
front ranks, crowding values, reference directions, and ideal/nadir points
directly.

See `API_MAPPING.md`, `ALGORITHM_NOTES.md`, and `VALIDATION.md` for details.
