# rgenoud-fortran

Modern Fortran/FPM translation of the computational core of the R package
`rgenoud` 5.9-0.3 (GENetic Optimization Using Derivatives).

The original package combines an evolutionary search with derivative-based
local refinement. This translation is standalone Fortran and has no R, Rcpp,
C, C++, BLAS, or LAPACK dependency.

## Main API

```fortran
use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize

type(genoud_options) :: opt
type(genoud_result) :: res
real(dp) :: lower(2), upper(2)

lower = -5.0_dp
upper =  5.0_dp
opt%pop_size = 100
opt%boundary_enforcement = 2

call genoud_optimize(objective, lower, upper, opt, res)
```

For vector-valued lexical optimization use `genoud_optimize_lexical`.

## Implemented computational features

- scalar minimization and maximization;
- the nine rgenoud operator weights, including cloning/survival accounting;
- uniform, boundary, and non-uniform mutation;
- polytope, simple, and heuristic crossover;
- whole non-uniform mutation;
- local-minimum crossover with configurable `P9` mixing;
- BFGS refinement of the current best individual;
- bounded projected-BFGS behavior for `boundary_enforcement = 2`;
- continuous and integer parameter modes;
- starting populations;
- MemoryMatrix-style exact duplicate evaluation cache;
- lexical/vector-valued objectives and custom lexical comparators;
- hard and soft generation limits;
- wait-generation and gradient convergence checks;
- analytical or numerical gradients;
- numerical Hessian output;
- population/sample moments.

## Build

```text
fpm build
fpm test
```

Examples:

```text
fpm run --example rosenbrock
fpm run --example integer_example
```

## Licensing and provenance

The upstream package declares `License: GPL-3`; this translation is therefore
released under GPL-3.0-only. The complete supplied upstream tree is retained in
`original/rgenoud-master/` for provenance and comparison.

See `TRANSLATION_COVERAGE.md` for deliberate runtime/API differences.

## Compiler portability

Version 0.1.2 routes user callbacks through typed module-level dispatch routines
whenever they would otherwise be invoked from nested internal procedures. This is
specifically compatible with gfortran/FPM builds that enable
`-Werror=implicit-interface`; it avoids relying on host-associated procedure
interfaces (optional or required) at nested call sites.
