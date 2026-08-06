! SPDX-License-Identifier: MIT
module uncorbets_types
  use uncorbets_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: uncorbets_ok = 0
  integer, parameter, public :: uncorbets_invalid_input = 1
  integer, parameter, public :: uncorbets_not_pos_semidefinite = 2
  integer, parameter, public :: uncorbets_singular_matrix = 3
  integer, parameter, public :: uncorbets_no_convergence = 4

  type, public :: status_type
    integer :: code = uncorbets_ok
    character(len=:), allocatable :: message
  contains
    procedure :: ok => status_ok
  end type status_type

  type, public :: torsion_result
    real(dp), allocatable :: matrix(:, :)
    integer :: iterations = 0
    real(dp) :: final_error = 0.0_dp
    logical :: converged = .true.
    type(status_type) :: status
  end type torsion_result

  type, public :: effective_bets_result
    real(dp), allocatable :: probability(:)
    real(dp) :: enb = 0.0_dp
    type(status_type) :: status
  end type effective_bets_result

  type, public :: max_effective_bets_result
    real(dp), allocatable :: weights(:)
    real(dp) :: enb = 0.0_dp
    integer :: objective_evaluations = 0
    integer :: gradient_evaluations = 0
    integer :: iterations = 0
    logical :: converged = .false.
    real(dp), allocatable :: lambda_lower(:)
    real(dp), allocatable :: lambda_upper(:)
    real(dp) :: lambda_equality = 0.0_dp
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: hessian(:, :)
    type(status_type) :: status
  end type max_effective_bets_result

  public :: set_status

contains

  logical function status_ok(self)
    class(status_type), intent(in) :: self
    status_ok = self%code == uncorbets_ok
  end function status_ok

  subroutine set_status(status, code, message)
    type(status_type), intent(out) :: status
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    status%code = code
    status%message = message
  end subroutine set_status

end module uncorbets_types
