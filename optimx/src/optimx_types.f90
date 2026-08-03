! SPDX-License-Identifier: GPL-2.0-only
module optimx_types
  use optimx_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: OPTIMX_SUCCESS = 0
  integer, parameter, public :: OPTIMX_MAXIT = 1
  integer, parameter, public :: OPTIMX_SMALL_GRADIENT = 2
  integer, parameter, public :: OPTIMX_BAD_UPDATE = 3
  integer, parameter, public :: OPTIMX_INVALID_INPUT = 20
  integer, parameter, public :: OPTIMX_BAD_EVALUATION = 21
  integer, parameter, public :: OPTIMX_LINESEARCH_FAILED = 22

  abstract interface
    subroutine optimx_callback(x, value, gradient, hessian, need_gradient, need_hessian, status)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      real(dp), intent(inout) :: gradient(:)
      real(dp), intent(inout) :: hessian(:, :)
      logical, intent(in) :: need_gradient, need_hessian
      integer, intent(out) :: status
    end subroutine optimx_callback
  end interface

  type, public :: optimx_problem
    integer :: n = 0
    procedure(optimx_callback), pointer, nopass :: objective => null()
    logical :: has_gradient = .false.
    logical :: has_hessian = .false.
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
    integer, allocatable :: mask(:)
  end type optimx_problem

  type, public :: optimx_control
    integer :: maxit = 500
    integer :: maxfeval = 5000
    integer :: trace = 0
    real(dp) :: reltol = 1.0e-8_dp
    real(dp) :: gradtol = 1.0e-7_dp
    real(dp) :: steptol = 1.0e-10_dp
    real(dp) :: acctol = 1.0e-4_dp
    real(dp) :: stepredn = 0.5_dp
    real(dp) :: initial_step = 1.0_dp
    logical :: maximize = .false.
    logical :: use_central = .true.
    logical :: kkt = .true.
    integer :: seed = 271828
  end type optimx_control

  type, public :: optimx_result
    real(dp), allocatable :: par(:)
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: hessian(:, :)
    real(dp) :: value = huge(1.0_dp)
    integer :: function_count = 0
    integer :: gradient_count = 0
    integer :: hessian_count = 0
    integer :: iterations = 0
    integer :: convergence = OPTIMX_INVALID_INPUT
    logical :: converged = .false.
    logical :: kkt1 = .false.
    logical :: kkt2 = .false.
    character(len=32) :: method = ''
    character(len=160) :: message = 'not run'
  end type optimx_result

  type, public :: optimx_multi_result
    type(optimx_result), allocatable :: runs(:)
    integer :: best = 0
  end type optimx_multi_result

  type, public :: derivative_check
    real(dp), allocatable :: analytic(:)
    real(dp), allocatable :: numeric(:)
    real(dp), allocatable :: error(:)
    real(dp) :: max_error = huge(1.0_dp)
    logical :: ok = .false.
    integer :: status = OPTIMX_INVALID_INPUT
  end type derivative_check

  type, public :: hessian_check
    real(dp), allocatable :: analytic(:, :)
    real(dp), allocatable :: numeric(:, :)
    real(dp), allocatable :: error(:, :)
    real(dp) :: max_error = huge(1.0_dp)
    logical :: ok = .false.
    integer :: status = OPTIMX_INVALID_INPUT
  end type hessian_check

  type, public :: kkt_result
    real(dp) :: projected_gradient_norm = huge(1.0_dp)
    real(dp) :: minimum_eigen_bound = -huge(1.0_dp)
    logical :: kkt1 = .false.
    logical :: kkt2 = .false.
    integer :: status = OPTIMX_INVALID_INPUT
  end type kkt_result

  type, public :: bounds_result
    real(dp), allocatable :: par(:)
    integer, allocatable :: mask(:)
    logical :: feasible = .false.
    integer :: status = OPTIMX_INVALID_INPUT
  end type bounds_result

  type, public :: scale_result
    real(dp) :: parameter_ratio = 1.0_dp
    real(dp) :: bound_ratio = 1.0_dp
    logical :: well_scaled = .true.
  end type scale_result

  type, public :: optsp_context
    real(dp) :: eps = sqrt(epsilon(1.0_dp))
    real(dp) :: fbase = 0.0_dp
  end type optsp_context

  public :: optimx_callback
end module optimx_types
