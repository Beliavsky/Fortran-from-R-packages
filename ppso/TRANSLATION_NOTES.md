# Translation notes

## Upstream

- Package: ppso
- Version: 0.9-99994
- DESCRIPTION date: 2024-11-20
- Author/Maintainer: Till Francke
- Language of computational source: R
- Declared license: Unlimited

## Source mapping

### `R/optim_pso.R` + `R/update_tasklist_pso.R`

Translated to `src/ppso_pso.f90`.

Preserved numerical details include:

- local and global best updates;
- asynchronous update-on-return behavior;
- synchronous whole-swarm behavior;
- velocity equation;
- the original whole-vector velocity limiter
  `V <- V * min(1, abs(Vmax/V))`, rather than componentwise clipping;
- parameter clipping after position updates;
- improvement/wait stopping criteria;
- initial estimates that remain fixed until first evaluation.

In synchronous mode, once an iteration begins the Fortran implementation
finishes the whole population, matching the R path even if this overshoots the
nominal function-call budget.

### `R/optim_dds.R` + `R/update_tasklist_dds.R`

Translated to `src/ppso_dds.f90`.

Preserved numerical details include:

- pre-search count `ceiling(max(0.005*max_calls, 5))`;
- decreasing DDS inclusion probability;
- at-least-one-dimension perturbation;
- normal perturbations with scale `r*(upper-lower)`;
- reflection followed by the package's over-reflection endpoint rule;
- all `part_xchange` modes.

The upstream serial implementation has a pre-search quirk: after separately
evaluating the first pre-search point, it removes that index before ranking and
before decrementing the function-call budget.  The Fortran control field
`legacy_serial_prerun_omission` reproduces the meaningful part of this behavior
when there are enough remaining pre-search points.  Corrected mode includes all
pre-search evaluations.

The R source can subsequently pair an omitted fitness with a different point's
coordinates in a narrow path because `X` is reset before the first
`update_tasklist_dds()` call while `fitness_X` is not.  This state corruption is
not reproduced.

### `R/init_particles.R`

Translated into `src/ppso_init.f90` and the optimizer state initialization.

- random box initialization is preserved;
- initial estimates use the first available particles and excess estimates are
  kept as pending evaluations;
- Latin-hypercube initialization is implemented natively, so the R `lhs`
  package is not required.

Project-file loading is represented by native Fortran checkpoint routines
rather than parsing R's TSV/data-frame serialization.

### Benchmark functions

`ackley_function.R`, `griewank_function.R`, `rastrigin_function.R`,
`sample_function.R`, and `sample_function2.R` are translated to
`src/ppso_benchmarks.f90`.  Artificial sleep/delay code used to test Rmpi is
omitted.

### Robust/Rmpi functions

`optim_ppso_robust.R`, `optim_pdds_robust.R`, `mpi_loop.R`,
`prepare_mpi_cluster.R`, `check_execution_timeout.R`, `close_mpi.R`, and
`request_object.R` primarily implement process management, retry/timeout
policy, message passing, and R object transport.  Those are not translated as
an MPI runtime.

Their underlying PSO/DDS update mathematics is shared with the serial routines
and is translated.  The state types are designed so a future MPI/coarray task
layer can evaluate requested points externally without changing those kernels.

### Plotting and visualization

`plot_optimization_progress.R`, `init_visualisation.R`, and
`do_plot_function.R` are omitted as requested.

## Fortran organization

- `ppso_kinds.f90`: floating-point kind
- `ppso_rng.f90`: deterministic RNG, normal draws, Latin hypercube
- `ppso_types.f90`: controls, results, optimizer states, callback interfaces
- `ppso_init.f90`: bounds and population initialization
- `ppso_pso.f90`: PSO engine
- `ppso_dds.f90`: DDS engine
- `ppso_checkpoint.f90`: state persistence
- `ppso_benchmarks.f90`: package objective examples
- `ppso.f90`: public facade module

All primary floating-point computations use

```fortran
integer, parameter :: dp = kind(1.0d0)
```

and all library procedures have explicit interfaces through modules.
