# API

## Modules

### `gradient`
Convenience umbrella module exporting the public optimizer API.

### `gradient_types`
Defines:

- `sqgde_options`
- `sqgde_result`
- `objective_fn`
- `get_algo_params`
- `validate_options`
- adaptation constants `SQGDE_RAND`, `SQGDE_CURRENT`, `SQGDE_BEST`
- convergence constants `CONVERGE_STDEV`, `CONVERGE_PERCENT`

### `gradient_sqgde`
Defines:

- `optim_sqgde`
- `adapt_sqgde_particle`
- `purify_population`

### `gradient_rng`
Defines `seed_rng` and internal random-number helpers.

## `sqgde_options`

Important members are:

```fortran
integer :: n_params
integer :: n_particles
integer :: n_diff
integer :: n_iter
real(dp), allocatable :: init_sd(:)
real(dp), allocatable :: init_center(:)
real(dp) :: step_size
real(dp) :: jitter_size
real(dp) :: crossover_rate
logical :: return_trace
integer :: thin
integer :: purify
integer :: adapt_scheme
integer :: give_up_init
integer :: stop_check
real(dp) :: stop_tol
integer :: converge_crit
```

`purify = 0` disables purification.  Positive values recompute every particle
objective on every `purify`-th iteration.

## `sqgde_result`

Contains:

```fortran
real(dp), allocatable :: solution(:)
real(dp) :: weight
logical :: converged
integer :: iterations
integer :: evaluations
integer :: status
character(len=160) :: message
real(dp), allocatable :: particles_trace(:,:,:)
real(dp), allocatable :: weights_trace(:,:)
integer :: trace_count
```

Traces are allocated only when `return_trace=.true.`.
