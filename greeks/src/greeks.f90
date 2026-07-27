! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
module greeks
  use greeks_kinds, only: dp
  use greeks_types
  use greeks_european, only: bs_european_greeks, bs_european_price
  use greeks_geometric_asian, only: bs_geometric_asian_greeks, bs_geometric_asian_price
  use greeks_american, only: binomial_american_greeks, binomial_american_price
  use greeks_malliavin, only: malliavin_european_greeks, malliavin_geometric_asian_greeks, &
    malliavin_asian_greeks, bs_malliavin_asian_greeks
  use greeks_implied_volatility, only: bs_implied_volatility, implied_volatility, &
    geometric_asian_implied_volatility, american_implied_volatility
  implicit none
  private
  public :: dp
  public :: greeks_result, mc_greeks_result, implied_vol_result
  public :: payoff_call, payoff_put, payoff_cash_call, payoff_cash_put
  public :: payoff_asset_call, payoff_asset_put
  public :: model_black_scholes, model_jump_diffusion
  public :: option_european, option_american, option_geometric_asian
  public :: payoff_callback, payoff_derivative_callback, jump_sampler_callback, price_callback
  public :: bs_european_greeks, bs_european_price
  public :: bs_geometric_asian_greeks, bs_geometric_asian_price
  public :: binomial_american_greeks, binomial_american_price
  public :: malliavin_european_greeks, malliavin_geometric_asian_greeks
  public :: malliavin_asian_greeks, bs_malliavin_asian_greeks
  public :: bs_implied_volatility, implied_volatility
  public :: geometric_asian_implied_volatility, american_implied_volatility
  public :: calculate_greeks
contains
  function calculate_greeks(spot, strike, rate, time, sigma, dividend, payoff, option_type, steps) result(res)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff, option_type
    integer, intent(in), optional :: steps
    type(greeks_result) :: res
    integer :: n
    select case (option_type)
    case (option_european)
      res=bs_european_greeks(spot,strike,rate,time,sigma,dividend,payoff)
    case (option_american)
      n=1000; if (present(steps)) n=steps
      res=binomial_american_greeks(spot,strike,rate,time,sigma,dividend,payoff,n)
    case (option_geometric_asian)
      res=bs_geometric_asian_greeks(spot,strike,rate,time,sigma,dividend,payoff)
    case default
      res%ok=.false.; res%message='unsupported option type'
    end select
  end function calculate_greeks
end module greeks
