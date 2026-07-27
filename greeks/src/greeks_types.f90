! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
module greeks_types
  use greeks_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: payoff_call = 1
  integer, parameter, public :: payoff_put = 2
  integer, parameter, public :: payoff_cash_call = 3
  integer, parameter, public :: payoff_cash_put = 4
  integer, parameter, public :: payoff_asset_call = 5
  integer, parameter, public :: payoff_asset_put = 6

  integer, parameter, public :: model_black_scholes = 1
  integer, parameter, public :: model_jump_diffusion = 2

  integer, parameter, public :: option_european = 1
  integer, parameter, public :: option_american = 2
  integer, parameter, public :: option_geometric_asian = 3

  type, public :: greeks_result
    real(dp) :: fair_value = 0.0_dp
    real(dp) :: delta = 0.0_dp
    real(dp) :: vega = 0.0_dp
    real(dp) :: theta = 0.0_dp
    real(dp) :: rho = 0.0_dp
    real(dp) :: epsilon = 0.0_dp
    real(dp) :: elasticity = 0.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: vanna = 0.0_dp
    real(dp) :: charm = 0.0_dp
    real(dp) :: vomma = 0.0_dp
    real(dp) :: veta = 0.0_dp
    real(dp) :: vera = 0.0_dp
    real(dp) :: speed = 0.0_dp
    real(dp) :: zomma = 0.0_dp
    real(dp) :: color = 0.0_dp
    real(dp) :: ultima = 0.0_dp
    real(dp) :: delta_d = 0.0_dp
    real(dp) :: rho_d = 0.0_dp
    real(dp) :: theta_d = 0.0_dp
    real(dp) :: vega_d = 0.0_dp
    real(dp) :: gamma_kombi = 0.0_dp
    logical :: ok = .true.
    character(len=160) :: message = ''
  end type greeks_result

  type, public :: mc_greeks_result
    type(greeks_result) :: estimate
    type(greeks_result) :: standard_error
    integer :: paths = 0
    integer :: steps = 0
    logical :: ok = .true.
    character(len=160) :: message = ''
  end type mc_greeks_result

  type, public :: implied_vol_result
    real(dp) :: volatility = 0.0_dp
    real(dp) :: price_error = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
    logical :: ok = .true.
    character(len=160) :: message = ''
  end type implied_vol_result

  abstract interface
    pure function payoff_callback(underlying, strike) result(value)
      import dp
      real(dp), intent(in) :: underlying, strike
      real(dp) :: value
    end function payoff_callback

    pure function payoff_derivative_callback(underlying, strike) result(value)
      import dp
      real(dp), intent(in) :: underlying, strike
      real(dp) :: value
    end function payoff_derivative_callback

    function jump_sampler_callback() result(value)
      import dp
      real(dp) :: value
    end function jump_sampler_callback

    function price_callback(volatility) result(value)
      import dp
      real(dp), intent(in) :: volatility
      real(dp) :: value
    end function price_callback
  end interface
  public :: payoff_callback, payoff_derivative_callback
  public :: jump_sampler_callback, price_callback
end module greeks_types
