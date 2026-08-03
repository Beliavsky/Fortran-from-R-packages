! SPDX-License-Identifier: LGPL-3.0-or-later
module nloptr_types
  use nloptr_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: NLOPT_FAILURE = -1
  integer, parameter, public :: NLOPT_INVALID_ARGS = -2
  integer, parameter, public :: NLOPT_OUT_OF_MEMORY = -3
  integer, parameter, public :: NLOPT_ROUNDOFF_LIMITED = -4
  integer, parameter, public :: NLOPT_FORCED_STOP = -5
  integer, parameter, public :: NLOPT_SUCCESS = 1
  integer, parameter, public :: NLOPT_STOPVAL_REACHED = 2
  integer, parameter, public :: NLOPT_FTOL_REACHED = 3
  integer, parameter, public :: NLOPT_XTOL_REACHED = 4
  integer, parameter, public :: NLOPT_MAXEVAL_REACHED = 5
  integer, parameter, public :: NLOPT_MAXTIME_REACHED = 6

  abstract interface
    subroutine objective_callback(x, value, gradient, need_gradient, status)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      real(dp), intent(inout) :: gradient(:)
      logical, intent(in) :: need_gradient
      integer, intent(out) :: status
    end subroutine objective_callback

    subroutine vector_callback(x, values, jacobian, need_jacobian, status)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: values(:)
      real(dp), intent(inout) :: jacobian(:, :)
      logical, intent(in) :: need_jacobian
      integer, intent(out) :: status
    end subroutine vector_callback
  end interface

  type, public :: nloptr_options
    character(len=40) :: algorithm = 'NLOPT_LD_LBFGS'
    character(len=40) :: local_algorithm = 'NLOPT_LD_LBFGS'
    real(dp) :: stopval = -huge(1.0_dp)
    real(dp) :: ftol_rel = 1.0e-8_dp
    real(dp) :: ftol_abs = 1.0e-10_dp
    real(dp) :: xtol_rel = 1.0e-7_dp
    real(dp) :: xtol_abs = 1.0e-9_dp
    real(dp) :: constraint_tol = 1.0e-7_dp
    real(dp) :: initial_step = 0.0_dp
    real(dp) :: penalty_initial = 10.0_dp
    real(dp) :: penalty_growth = 10.0_dp
    real(dp) :: max_time = 0.0_dp
    integer :: maxeval = 2000
    integer :: max_outer = 8
    integer :: population = 0
    integer :: seed = 271828
    integer :: print_level = 0
    logical :: check_derivatives = .false.
  end type nloptr_options

  type, public :: nloptr_problem
    integer :: n = 0
    integer :: n_ineq = 0
    integer :: n_eq = 0
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
    procedure(objective_callback), pointer, nopass :: objective => null()
    procedure(vector_callback), pointer, nopass :: inequality => null()
    procedure(vector_callback), pointer, nopass :: equality => null()
  end type nloptr_problem

  type, public :: nloptr_result
    real(dp), allocatable :: solution(:)
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: max_constraint = huge(1.0_dp)
    integer :: status = NLOPT_FAILURE
    integer :: iterations = 0
    integer :: evaluations = 0
    logical :: converged = .false.
    character(len=160) :: message = 'not run'
  end type nloptr_result

  type, public :: derivative_check_result
    real(dp), allocatable :: analytic(:, :)
    real(dp), allocatable :: numeric(:, :)
    real(dp), allocatable :: relative_error(:, :)
    logical, allocatable :: warning(:, :)
    integer :: n_warnings = 0
    integer :: status = 0
  end type derivative_check_result

  public :: objective_callback, vector_callback
end module nloptr_types
