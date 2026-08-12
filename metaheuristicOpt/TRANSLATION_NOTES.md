# Translation notes

## Scope

This release translates the computational content of `metaheuristicOpt` 2.0.0
from R to modern Fortran.  The upstream package has no C/C++/Fortran backend;
all 21 optimization engines are implemented in R files under `R/`.

Translated components:

- population generation and bound clipping;
- objective transformation for minimization/maximization;
- roulette selection;
- all 21 optimization engines;
- algorithm-specific random walks, Levy flights, selection/crossover/mutation,
  social-force calculations, scout/migration-like replacement, and histories;
- `metaOpt`-style name dispatch;
- benchmark/helper functions in `metaheuristic.FunctionCollection.R`.

Not translated because they are R/UI infrastructure rather than numerical
algorithms:

- `txtProgressBar` and console progress reporting;
- commented plotting calls and graphics imports;
- R list/data-frame/matrix presentation conventions;
- roxygen/man-page machinery and timing returned by `system.time()`.

## Data model

R matrices whose rows are candidate solutions become `real(dp)` arrays with
shape `(population, variables)`.  Bounds are explicit vectors.  The objective is
a Fortran callback with interface

```fortran
real(dp) function f(x)
   real(dp), intent(in) :: x(:)
end function
```

A private Park-Miller generator supplies deterministic uniform random numbers;
normal values use Box-Muller where Levy-flight code requires them.  Consequently
an integer seed is reproducible within this Fortran project but is not expected
to reproduce R's Mersenne-Twister stream.

## Fidelity choices

The goal is source-level behavioral fidelity, including unusual formulas where
they are part of the active package code.  Examples include:

- PSO objective evaluation after each coordinate update;
- GA roulette weights based on the package's internal fitness values;
- GOA's pairwise two-coordinate social interaction and dummy coordinate for odd
  dimensions;
- DA's componentwise neighbor test and package Levy construction;
- CSO's seeking-copy multiplication expression;
- KH's reciprocal-fitness food location and genetic operators;
- ABC's initialization over `min(lower):max(upper)`;
- DE's mutually exclusive permutation vectors and the source spelling
  `"clasical"`.

`mh_control%legacy_quirks` isolates source behaviors that are useful for
compatibility but are almost certainly mistakes.  Corrected behavior is
available by setting it false.  See README.md for the complete list.

## Known deviations

A few R implementation accidents are intentionally normalized where literal
reproduction would be unsafe or not meaningful in Fortran:

- invalid/non-integer indexing is replaced by explicit integer sampling;
- impossible population sizes produce `error stop` rather than potentially
  entering R rejection loops forever;
- dimensions and bounds are validated explicitly;
- the Fortran result always reports an explicit objective-evaluation count.

The R `metaOpt()` wrapper sometimes fails to forward algorithm-specific values
from its `control` list and instead lets the direct algorithm wrapper use its
defaults.  The Fortran `metaopt()` dispatcher passes the corresponding
`mh_control` fields deliberately; direct algorithm routines and dispatcher calls
therefore share one consistent control model.

## Files

- `src/mh_support.f90`: kind, RNG, shuffling, sampling, roulette utilities.
- `src/metaheuristic_opt.f90`: public types, dispatcher, helpers, and 21 engines.
- `test/`: regression/smoke tests.
- `example/`: FPM examples.
- `original/metaheuristicOpt-master/`: unmodified attached package source.
