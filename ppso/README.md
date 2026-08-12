# ppso-fortran 0.1.0

Modern Fortran translation of the computational core of the R package **ppso
0.9-99994** by Till Francke.

The original package implements Particle Swarm Optimization (PSO) and
Dynamically Dimensioned Search (DDS), with optional Rmpi-based parallel task
execution and project-file restart support.  This translation keeps the
numerical search algorithms and state required for restart, but does not
reimplement R, Rmpi, plotting, or GUI/interactive infrastructure.

## Implemented computational functionality

- Particle Swarm Optimization (`optim_pso`)
  - asynchronous serial update mode, matching `wait_complete_iteration=FALSE`
  - synchronous whole-swarm mode, matching `wait_complete_iteration=TRUE`
  - inertia (`w`), cognitive (`c1`), and social (`c2`) terms
  - ppso's whole-vector `Vmax` scaling rule
  - bound clipping
  - local/global best state
  - absolute/relative improvement and wait stopping rules
  - maximum-iteration and function-call stopping
  - supplied initial estimates, including pending estimates
- Dynamically Dimensioned Search (`optim_dds`)
  - Tolson/Shoemaker decreasing dimension-inclusion probability
  - Gaussian perturbations scaled by parameter range and `r`
  - ppso's reflection/over-reflection boundary rule
  - multiple DDS particles
  - all four `part_xchange=0..3` communication/relocation modes
  - ppso pre-search initialization
  - optional compatibility with the serial R pre-search omission quirk
- Random and native Latin-hypercube initialization
- Checkpoint serialization for PSO and DDS state
- Original package benchmark functions:
  - `rastrigin_function`
  - `ackley_function`
  - `griewank_function`
  - `sample_function`
  - `sample_function2`
- Optional monitor callback
- Deterministic standalone RNG

## Not translated

The following are runtime/interface infrastructure rather than the numerical
optimization algorithms and are intentionally omitted:

- Rmpi process creation, message passing, worker timeout handling, and remote
  working-directory management
- `rgl` and base-R optimization visualizations
- interactive `readline()` debugging
- R list/data-frame/formula machinery
- R logfile formatting and the exact tab-separated R project-file format

The asynchronous serial PSO path still models the important algorithmic effect
of result order: after one particle returns, its result can update the global
best before the next particle is evaluated.  Native parallel execution can be
added around the state API later without changing the optimization kernels.

## Build with FPM

```text
fpm build
fpm test
fpm run --example pso_rastrigin
fpm run --example dds_ackley
```

The package has no external numerical-library dependency.

FPM was not installed in the translation environment, so the FPM layout was
validated using equivalent direct GNU Fortran 14.2.0 builds.

## Basic PSO example

```fortran
program example
    use ppso, only : dp, pso_control, ppso_result, optim_pso, rastrigin_function
    implicit none

    type(pso_control) :: control
    type(ppso_result) :: result
    real(dp) :: lower(2), upper(2)

    lower = -1.0_dp
    upper =  1.0_dp

    control%number_of_particles = 40
    control%max_number_of_iterations = 80
    control%max_number_function_calls = 3200
    control%seed = 20260810_8

    call optim_pso(rastrigin_function, lower, upper, result, control)
    print *, result%value
    print *, result%par
end program example
```

The ppso Rastrigin function is

```text
sum(x**2 - cos(2*pi*x))
```

so its two-dimensional global minimum is `-2` at `(0,0)`.

## DDS example

```fortran
use ppso, only : dp, dds_control, ppso_result, optim_dds, ackley_function

type(dds_control) :: control
type(ppso_result) :: result
real(dp) :: lower(4), upper(4)

lower = -5.0_dp
upper =  5.0_dp
control%number_of_particles = 2
control%max_number_function_calls = 1600
control%part_xchange = 2
control%legacy_serial_prerun_omission = .false.

call optim_dds(ackley_function, lower, upper, result, control)
```

## DDS serial pre-search compatibility

The upstream serial `optim_dds.R` evaluates the first pre-search point
separately, then removes that point from the vector used for ranking and call
accounting.  For the common one-particle case this means one real objective
call is omitted from the reported count and cannot become the selected
pre-search best.

For source compatibility the Fortran default is:

```fortran
control%legacy_serial_prerun_omission = .true.
```

Set it to `.false.` for the corrected behavior.  `ppso_result` contains both:

- `function_calls`: package-compatible reported count
- `actual_function_calls`: actual objective callback count

The translation does **not** reproduce a more severe downstream R state
inconsistency that can associate the omitted point's fitness with another
point's coordinates.

## Checkpointing

The public routines

```text
save_pso_state / load_pso_state
save_dds_state / load_dds_state
```

serialize the complete numerical state, including RNG state, so a Fortran
application can persist and resume a search reproducibly.  The binary format is
specific to this Fortran port; it is not the original R tab-separated project
file format.

## Random-number compatibility

The translation uses a self-contained xorshift RNG plus Box-Muller normals.
The same integer seed is reproducible within this Fortran implementation, but
it does not generate R's random stream, so trajectories are not expected to be
bitwise identical to R.

## Validation

The regression suite covers:

1. all original benchmark-function optima at zero;
2. asynchronous and synchronous PSO on the package Rastrigin function;
3. DDS on a smooth sphere objective using all `part_xchange` modes;
4. legacy and corrected serial DDS pre-search call accounting;
5. deterministic checkpoint save/load equivalence.

The suite is tested with both optimized compilation and a debug build using
bounds checking, warnings, and implicit-interface errors.

## License and provenance

The upstream `DESCRIPTION` declares:

```text
License: Unlimited
```

No separate license text is present in the supplied repository.  This
translation preserves that declaration and retains the complete supplied R
package under `original/` for attribution and auditability.  See `LICENSES.md`
and `TRANSLATION_NOTES.md`.
