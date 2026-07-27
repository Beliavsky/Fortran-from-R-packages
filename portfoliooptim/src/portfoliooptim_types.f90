! SPDX-License-Identifier: GPL-3.0-only
! Based on PortfolioOptim 1.1.1 by Andrzej Palczewski and Aleksandra Dabrowska.
module portfoliooptim_types
  use portfoliooptim_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: risk_cvar = 1
  integer, parameter, public :: risk_dcvar = 2
  integer, parameter, public :: risk_lsad = 3
  integer, parameter, public :: risk_mad = 4

  type, public :: risk_result
    real(dp) :: var = 0.0_dp
    real(dp) :: cvar = 0.0_dp
    real(dp) :: mad = 0.0_dp
    real(dp) :: mean = 0.0_dp
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type risk_result

  type, public :: lp_result
    real(dp), allocatable :: x(:)
    real(dp) :: objective = huge(1.0_dp)
    integer :: iterations = 0
    logical :: optimal = .false.
    logical :: feasible = .false.
    logical :: unbounded = .false.
    character(len=160) :: message = ''
  end type lp_result

  type, public :: projection_result
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
    real(dp) :: residual = huge(1.0_dp)
    integer :: iterations = 0
    logical :: converged = .false.
    character(len=160) :: message = ''
  end type projection_result

  type, public :: portfolio_result
    real(dp), allocatable :: return_mean(:)
    real(dp), allocatable :: theta(:)
    real(dp) :: mu = 0.0_dp
    real(dp) :: cvar = 0.0_dp
    real(dp) :: var = 0.0_dp
    real(dp) :: mad = 0.0_dp
    real(dp) :: risk = 0.0_dp
    real(dp) :: new_portfolio_return = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
    character(len=160) :: message = ''
  end type portfolio_result

end module portfoliooptim_types
