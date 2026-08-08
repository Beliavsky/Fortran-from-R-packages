! SPDX-License-Identifier: GPL-3.0-only
module bb_types
  use bb_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: bb_success = 0
  integer, parameter, public :: bb_max_iterations = 1
  integer, parameter, public :: bb_max_evaluations = 2
  integer, parameter, public :: bb_function_error = 3
  integer, parameter, public :: bb_gradient_error = 4
  integer, parameter, public :: bb_projection_error = 5
  integer, parameter, public :: bb_stagnation = 6
  integer, parameter, public :: bb_line_search_error = 7
  integer, parameter, public :: bb_anomalous_iteration = 8
  integer, parameter, public :: bb_no_improvement = 9
  integer, parameter, public :: bb_invalid_input = 10

  type, public :: spg_control
    integer :: m = 10
    integer :: maxit = 1500
    integer :: maxfeval = 10000
    integer :: method = 3
    real(dp) :: ftol = 1.0e-10_dp
    real(dp) :: gtol = 1.0e-5_dp
    real(dp) :: eps = 1.0e-7_dp
    logical :: maximize = .false.
    logical :: trace = .false.
    integer :: triter = 10
  end type spg_control

  type, public :: spg_result
    real(dp), allocatable :: par(:)
    real(dp) :: value = huge(1.0_dp)
    real(dp) :: gradient = huge(1.0_dp)
    real(dp) :: fn_reduction = 0.0_dp
    integer :: iter = 0
    integer :: feval = 0
    integer :: convergence = bb_invalid_input
    integer :: method = 0
    integer :: m = 0
    character(len=:), allocatable :: message
  contains
    procedure :: succeeded => spg_succeeded
  end type spg_result

  type, public :: sane_control
    integer :: maxit = 1500
    integer :: m = 10
    integer :: method = 2
    real(dp) :: tol = 1.0e-7_dp
    logical :: trace = .false.
    integer :: triter = 10
    integer :: noimp = 100
    logical :: nm = .false.
    logical :: bfgs = .false.
  end type sane_control

  type, public :: sane_result
    real(dp), allocatable :: par(:)
    real(dp) :: residual = huge(1.0_dp)
    real(dp) :: fn_reduction = 0.0_dp
    integer :: feval = 0
    integer :: iter = 0
    integer :: convergence = bb_invalid_input
    integer :: method = 0
    integer :: m = 0
    logical :: nm_used = .false.
    character(len=:), allocatable :: message
  contains
    procedure :: succeeded => sane_succeeded
  end type sane_result

  type, public :: bboptim_control
    integer :: maxit = 1500
    integer :: methods(3) = [2, 3, 1]
    integer :: m_values(2) = [50, 10]
    real(dp) :: ftol = 1.0e-10_dp
    real(dp) :: gtol = 1.0e-5_dp
    integer :: maxfeval = 10000
    logical :: maximize = .false.
    logical :: trace = .false.
    integer :: triter = 10
    real(dp) :: eps = 1.0e-7_dp
  end type bboptim_control

  type, public :: bbsolve_control
    integer :: maxit = 1500
    integer :: methods(3) = [2, 3, 1]
    integer :: m_values(2) = [50, 10]
    real(dp) :: tol = 1.0e-7_dp
    logical :: trace = .false.
    integer :: triter = 10
    integer :: noimp = 100
    logical :: try_nm = .true.
  end type bbsolve_control

  type, public :: multistart_result
    real(dp), allocatable :: par(:, :)
    real(dp), allocatable :: fvalue(:)
    logical, allocatable :: converged(:)
  end type multistart_result

contains

  logical function spg_succeeded(self)
    class(spg_result), intent(in) :: self
    spg_succeeded = self%convergence == bb_success
  end function spg_succeeded

  logical function sane_succeeded(self)
    class(sane_result), intent(in) :: self
    sane_succeeded = self%convergence == bb_success
  end function sane_succeeded
end module bb_types
