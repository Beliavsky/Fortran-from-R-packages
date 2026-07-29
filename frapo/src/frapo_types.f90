! SPDX-License-Identifier: GPL-3.0-or-later
module frapo_types
  use frapo_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: frapo_ok = 0
  integer, parameter, public :: frapo_invalid_input = 1
  integer, parameter, public :: frapo_singular = 2
  integer, parameter, public :: frapo_no_convergence = 3
  integer, parameter, public :: frapo_infeasible = 4

  type, public :: portfolio_result
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: drawdowns(:)
    real(dp) :: objective = 0.0_dp
    real(dp) :: terminal_return = 0.0_dp
    real(dp) :: threshold = 0.0_dp
    real(dp) :: risk_value = 0.0_dp
    integer :: status = frapo_ok
    integer :: iterations = 0
    character(len=64) :: portfolio_type = ''
  end type portfolio_result

  type, public :: optimizer_result
    real(dp), allocatable :: x(:)
    real(dp) :: objective = 0.0_dp
    real(dp) :: primal_residual = huge(1.0_dp)
    real(dp) :: dual_residual = huge(1.0_dp)
    integer :: status = frapo_ok
    integer :: iterations = 0
  end type optimizer_result
end module frapo_types
