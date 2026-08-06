! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

module rpeglmen_types
  use rpeglmen_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: rpe_success = 0
  integer, parameter, public :: rpe_invalid_input = 1
  integer, parameter, public :: rpe_no_convergence = 2
  integer, parameter, public :: rpe_line_search_failure = 3
  integer, parameter, public :: rpe_numerical_failure = 4

  integer, parameter, public :: model_exponential = 1
  integer, parameter, public :: model_gamma = 2

  type, public :: enet_options
    real(dp) :: alpha = 0.5_dp
    integer :: num_lambda = 100
    real(dp) :: min_lambda_ratio = 1.0e-4_dp
    real(dp) :: min_lambda_absolute = 0.0_dp
    integer :: max_iter = 1000
    real(dp) :: abs_tol = 1.0e-8_dp
    real(dp) :: rel_tol = 1.0e-8_dp
    real(dp) :: initial_step = 1.0_dp
    real(dp) :: backtrack = 0.5_dp
    integer :: max_backtrack = 80
    integer :: k_fold = 5
    integer :: k_fold_iter = 5
    integer :: seed = 12345
    logical :: has_intercept = .true.
    logical :: penalize_intercept = .false.
    logical :: normalize_gradient = .false.
    logical :: use_fista = .true.
    logical :: source_proximal = .false.
    character(len=16) :: cv_metric = 'nll'
  end type enet_options

  type, public :: fit_result
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: objective(:)
    real(dp), allocatable :: lambda_grid(:)
    real(dp), allocatable :: cv_mean(:)
    real(dp), allocatable :: cv_sd(:)
    real(dp) :: shape = 1.0_dp
    real(dp) :: selected_lambda = 0.0_dp
    integer :: iterations = 0
    integer :: status = rpe_success
    logical :: converged = .false.
    character(len=160) :: message = ''
  end type fit_result

  public :: clear_fit_result

  type, public :: path_result
    real(dp), allocatable :: coefficients(:, :)
    real(dp), allocatable :: lambda_grid(:)
    integer, allocatable :: iterations(:)
    logical, allocatable :: converged(:)
    integer :: status = rpe_success
    character(len=160) :: message = ''
  end type path_result

contains

  subroutine clear_fit_result(result)
    type(fit_result), intent(inout) :: result

    if (allocated(result%coefficients)) deallocate(result%coefficients)
    if (allocated(result%objective)) deallocate(result%objective)
    if (allocated(result%lambda_grid)) deallocate(result%lambda_grid)
    if (allocated(result%cv_mean)) deallocate(result%cv_mean)
    if (allocated(result%cv_sd)) deallocate(result%cv_sd)
    result%shape = 1.0_dp
    result%selected_lambda = 0.0_dp
    result%iterations = 0
    result%status = rpe_success
    result%converged = .false.
    result%message = ''
  end subroutine clear_fit_result

end module rpeglmen_types
