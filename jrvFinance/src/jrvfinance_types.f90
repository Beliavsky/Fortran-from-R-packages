! SPDX-License-Identifier: GPL-2.0-or-later
module jrvfinance_types
  use jrvfinance_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: JRV_OK = 0
  integer, parameter, public :: JRV_INVALID_ARGUMENT = 1
  integer, parameter, public :: JRV_NO_CONVERGENCE = 2
  integer, parameter, public :: JRV_NO_ROOT = 3
  integer, parameter, public :: JRV_NONFINITE = 4
  integer, parameter, public :: JRV_SIZE_MISMATCH = 5
  real(dp), parameter, public :: CONTINUOUS_FREQUENCY = 0.0_dp

  type, public :: date_t
    integer :: year = 1970
    integer :: month = 1
    integer :: day = 1
  end type date_t

  type, public :: root_result
    real(dp) :: root = 0.0_dp
    real(dp) :: value = 0.0_dp
    integer :: iterations = 0
    integer :: status = JRV_OK
    character(len=160) :: message = ''
  end type root_result

  type, public :: annuity_breakup_result
    real(dp) :: opening_principal = 0.0_dp
    real(dp) :: interest_part = 0.0_dp
    real(dp) :: principal_part = 0.0_dp
    real(dp) :: closing_principal = 0.0_dp
    integer :: status = JRV_OK
  end type annuity_breakup_result

  type, public :: bond_cashflows
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: cashflow(:)
    real(dp) :: accrued = 0.0_dp
    integer :: status = JRV_OK
  end type bond_cashflows

  type, public :: black_scholes_result
    real(dp) :: call = 0.0_dp
    real(dp) :: put = 0.0_dp
    real(dp) :: call_delta = 0.0_dp
    real(dp) :: put_delta = 0.0_dp
    real(dp) :: call_theta = 0.0_dp
    real(dp) :: put_theta = 0.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: vega = 0.0_dp
    real(dp) :: call_rho = 0.0_dp
    real(dp) :: put_rho = 0.0_dp
    real(dp) :: d1 = 0.0_dp
    real(dp) :: d2 = 0.0_dp
    real(dp) :: nd1 = 0.0_dp
    real(dp) :: nd2 = 0.0_dp
    real(dp) :: nminusd1 = 0.0_dp
    real(dp) :: nminusd2 = 0.0_dp
    real(dp) :: call_probability = 0.0_dp
    real(dp) :: put_probability = 0.0_dp
    integer :: status = JRV_OK
  end type black_scholes_result

end module jrvfinance_types
