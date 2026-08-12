# GA-fortran

Modern Fortran translation of the computational algorithms in the R package
**GA 3.2.5** by Luca Scrucca.

The original package is licensed under GPL (>= 2). This translation is
therefore distributed under **GPL-2.0-or-later**. See `LICENSE`, `LICENSES.md`,
and `original/`.

## Scope

The translation concentrates on numerical and stochastic optimization code.
R S4 classes, plotting, printing, console monitors, `foreach`/cluster setup,
Rcpp marshalling, and other R-runtime infrastructure are intentionally omitted.

Implemented package algorithms include:

- generational genetic algorithms for
  - binary chromosomes,
  - bounded real-valued chromosomes,
  - permutation chromosomes;
- suggestions/seed populations;
- fitness caching and evaluation counts;
- GA fitness summaries and `run` / `maxFitness` stopping;
- elitism;
- adaptive mutation probability (`ga_pmutation`);
- Differential Evolution (`de` / `gareal_de` semantics);
- coarse-grained island GAs with ring migration;
- binary/decimal and binary/Gray conversions;
- bound repair and reflection helpers;
- all native selection, crossover, and mutation families listed below.

## Genetic operators

### Selection

The constants exposed by module `ga` are:

- `SEL_LINEAR_RANK` -- linear-rank selection (`ga_lrSelection`)
- `SEL_NONLINEAR_RANK` -- nonlinear-rank selection (`ga_nlrSelection`)
- `SEL_ROULETTE` -- roulette-wheel selection (`ga_rwSelection`)
- `SEL_TOURNAMENT` -- tournament selection (`ga_tourSelection`)
- `SEL_REAL_LINEAR_SCALE` -- real-valued linear fitness scaling
  (`gareal_lsSelection`), the GA package default for real chromosomes
- `SEL_REAL_SIGMA` -- sigma scaling (`gareal_sigmaSelection`)

### Crossover

- `CROSS_SINGLE_POINT` -- generic single-point crossover
- `CROSS_BINARY_UNIFORM` -- binary uniform crossover
- `CROSS_REAL_WEIGHTED` -- whole arithmetic crossover
- `CROSS_REAL_LOCAL` -- local arithmetic crossover, the default real crossover
- `CROSS_REAL_BLX` -- BLX-alpha crossover
- `CROSS_REAL_LAPLACE` -- Laplace crossover
- `CROSS_PERM_CYCLE` -- cycle crossover
- `CROSS_PERM_PMX` -- partially matched crossover
- `CROSS_PERM_OX` -- order crossover, the default permutation crossover
- `CROSS_PERM_PBX` -- position-based crossover

### Mutation

- `MUT_BINARY_RANDOM` -- random bit flip
- `MUT_REAL_RANDOM` -- uniform random mutation, the default real mutation
- `MUT_REAL_NONUNIFORM` -- non-uniform mutation
- `MUT_REAL_RANDOM_SHIFT` -- random-shift mutation
- `MUT_REAL_POWER` -- power mutation
- `MUT_PERM_INVERSION` -- simple inversion, the default permutation mutation
- `MUT_PERM_INSERTION` -- insertion mutation
- `MUT_PERM_SWAP` -- swap mutation
- `MUT_PERM_DISPLACEMENT` -- displacement mutation
- `MUT_PERM_SCRAMBLE` -- scramble mutation

## Basic real-valued GA

```fortran
program demo
  use ga
  implicit none
  type(ga_control_type) :: control
  type(ga_real_result) :: result
  real(dp) :: lower(3), upper(3)

  lower = -5.0_dp
  upper =  5.0_dp
  control%pop_size = 80
  control%max_iter = 300
  control%seed = 1234

  call ga_real(fitness, lower, upper, control, result)
  print *, result%fitness_value
  print *, result%solution
contains
  function fitness(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = -sum((x-[1.0_dp,-2.0_dp,0.5_dp])**2)
  end function fitness
end program demo
```

`ga_binary` and `ga_permutation` use integer-vector fitness callbacks. A generic
interface `ga_optimize` dispatches to the three routines from the argument
shapes.

## Differential Evolution

`de_real` follows the DE operator in `R/gade.R` and `src/genope.cpp`:

```fortran
type(de_control_type) :: control
control%stepsize = 0.8_dp
control%pcrossover = 0.5_dp
call de_real(fitness, lower, upper, control, result)
```

Set `control%dither=.true.` to draw the scale factor independently in `[0.5,1]`
for each target vector, corresponding to the package behavior when `F` is `NA`.

## Islands

The serial routines

- `ga_islands_real`
- `ga_islands_binary`
- `ga_islands_permutation`

implement the package's coarse-grained ring migration. The best fraction of
one island migrates to the next and replaces random non-elite individuals.
The original R package can distribute islands using `foreach`; parallel
runtime orchestration is intentionally left to the Fortran application.

## R features intentionally not translated

The following are interface/orchestration layers rather than package-owned
numerical kernels:

- S4 `ga`, `de`, and `gaisl` classes and methods;
- plotting, palettes, perspective plots, printing, and monitors;
- R formula/list argument handling;
- `foreach`, `doParallel`, and `doRNG` cluster management;
- calls to `stats::optim` for optional hybrid local search;
- `updatePop`'s R-attribute return convention and arbitrary `postFitness`
  mutation of S4 objects.

The GA likelihood/fitness callback is native Fortran, so an application can
call any desired Fortran local optimizer around `ga_real` if hybrid search is
needed.

## Reproducibility

`control%seed` initializes Fortran's intrinsic RNG deterministically in this
implementation. It does **not** reproduce R's RNG stream for the same numeric
seed. Algorithmic distributions and operator semantics are preserved, not
cross-language bit-identical random draws.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example real_ga_example
```

No external numerical libraries are required.

The translation was also validated directly with GNU Fortran 14.2.0 using
Fortran 2018, runtime bounds checking, warnings, and optimized builds.

## Source layout

- `src/ga_base.f90` -- kinds, RNG, utility/statistical helpers
- `src/ga_operators.f90` -- GA populations and genetic operators
- `src/ga_core.f90` -- binary/real/permutation GA and DE engines
- `src/ga_islands.f90` -- coarse-grained islands
- `src/ga.f90` -- public convenience module
- `test/` -- regression and invariant tests
- `example/` -- small executable examples
- `original/` -- original package metadata and computational sources retained
  for license/provenance auditing

See `TRANSLATION_NOTES.md` for a more detailed R-to-Fortran mapping.
