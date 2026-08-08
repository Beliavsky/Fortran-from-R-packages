! SPDX-License-Identifier: MIT
module gradient_types
  use gradient_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: SQGDE_RAND = 1
  integer, parameter, public :: SQGDE_CURRENT = 2
  integer, parameter, public :: SQGDE_BEST = 3
  integer, parameter, public :: CONVERGE_STDEV = 1
  integer, parameter, public :: CONVERGE_PERCENT = 2

  type, public :: sqgde_options
    integer :: n_params = 0
    integer :: n_particles = 0
    integer :: n_diff = 2
    integer :: n_iter = 1000
    real(dp), allocatable :: init_sd(:)
    real(dp), allocatable :: init_center(:)
    real(dp) :: step_size = -1.0_dp
    real(dp) :: jitter_size = 1.0e-6_dp
    real(dp) :: crossover_rate = 1.0_dp
    logical :: return_trace = .false.
    integer :: thin = 1
    integer :: purify = 0
    integer :: adapt_scheme = SQGDE_RAND
    integer :: give_up_init = 100
    integer :: stop_check = 10
    real(dp) :: stop_tol = 1.0e-4_dp
    integer :: converge_crit = CONVERGE_STDEV
  end type sqgde_options

  type, public :: sqgde_result
    real(dp), allocatable :: solution(:)
    real(dp) :: weight = huge(1.0_dp)
    logical :: converged = .false.
    integer :: iterations = 0
    integer :: evaluations = 0
    integer :: status = 0
    character(len=160) :: message = ''
    real(dp), allocatable :: particles_trace(:,:,:)
    real(dp), allocatable :: weights_trace(:,:)
    integer :: trace_count = 0
  end type sqgde_result

  abstract interface
    function objective_fn(x) result(f)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function objective_fn
  end interface
  public :: objective_fn

  public :: get_algo_params, validate_options
contains
  subroutine get_algo_params(n_params, options)
    integer, intent(in) :: n_params
    type(sqgde_options), intent(out) :: options

    options%n_params = n_params
    options%n_particles = max(3*n_params, 4)
    options%n_diff = 2
    options%n_iter = 1000
    allocate(options%init_sd(n_params), options%init_center(n_params))
    options%init_sd = 0.01_dp
    options%init_center = 0.0_dp
    options%step_size = 2.38_dp/sqrt(2.0_dp*real(max(n_params,1),dp))
    options%jitter_size = 1.0e-6_dp
    options%crossover_rate = 1.0_dp
    options%return_trace = .false.
    options%thin = 1
    options%purify = 0
    options%adapt_scheme = SQGDE_RAND
    options%give_up_init = 100
    options%stop_check = 10
    options%stop_tol = 1.0e-4_dp
    options%converge_crit = CONVERGE_STDEV
  end subroutine get_algo_params

  subroutine validate_options(options, status, message)
    type(sqgde_options), intent(in) :: options
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    status = 0
    message = ''
    if (options%n_params < 1) then
      status = 1; message = 'n_params must be positive'; return
    end if
    if (options%n_particles < 4) then
      status = 2; message = 'n_particles must be at least 4'; return
    end if
    if (options%n_iter < 4) then
      status = 3; message = 'n_iter must be at least 4'; return
    end if
    if (.not. allocated(options%init_sd) .or. size(options%init_sd) /= options%n_params) then
      status = 4; message = 'init_sd must have n_params elements'; return
    end if
    if (.not. allocated(options%init_center) .or. size(options%init_center) /= options%n_params) then
      status = 5; message = 'init_center must have n_params elements'; return
    end if
    if (any(options%init_sd <= 0.0_dp)) then
      status = 6; message = 'init_sd must be positive'; return
    end if
    if (options%step_size <= 0.0_dp) then
      status = 7; message = 'step_size must be positive'; return
    end if
    if (options%jitter_size < 0.0_dp) then
      status = 8; message = 'jitter_size cannot be negative'; return
    end if
    if (options%crossover_rate <= 0.0_dp .or. options%crossover_rate > 1.0_dp) then
      status = 9; message = 'crossover_rate must be in (0,1]'; return
    end if
    if (options%n_diff < 1) then
      status = 10; message = 'n_diff must be positive'; return
    end if
    if (2*options%n_diff + 1 > options%n_particles) then
      status = 11
      message = 'SQG-DE requires 2*n_diff+1 <= n_particles for unique parent sampling'
      return
    end if
    if (options%thin < 1) then
      status = 12; message = 'thin must be positive'; return
    end if
    if (options%purify < 0) then
      status = 13; message = 'purify must be zero (disabled) or positive'; return
    end if
    if (options%adapt_scheme < SQGDE_RAND .or. options%adapt_scheme > SQGDE_BEST) then
      status = 14; message = 'invalid adapt_scheme'; return
    end if
    if (options%give_up_init < 1) then
      status = 15; message = 'give_up_init must be positive'; return
    end if
    if (options%stop_check < 2) then
      status = 16; message = 'stop_check must be at least 2'; return
    end if
    if (options%stop_tol < 0.0_dp) then
      status = 17; message = 'stop_tol cannot be negative'; return
    end if
    if (options%converge_crit < CONVERGE_STDEV .or. options%converge_crit > CONVERGE_PERCENT) then
      status = 18; message = 'invalid convergence criterion'; return
    end if
  end subroutine validate_options
end module gradient_types
