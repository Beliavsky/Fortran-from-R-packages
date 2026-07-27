! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
module greeks_american
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use greeks_kinds, only: dp
  use greeks_types, only: greeks_result, payoff_call, payoff_put
  use greeks_european, only: bs_european_price
  use greeks_math, only: safe_divide
  implicit none
  private
  public :: binomial_american_greeks, binomial_american_price
contains
  pure elemental function intrinsic_value(spot, strike, payoff) result(value)
    real(dp), intent(in) :: spot, strike
    integer, intent(in) :: payoff
    real(dp) :: value
    if (payoff == payoff_call) then
      value = max(spot-strike, 0.0_dp)
    else
      value = max(strike-spot, 0.0_dp)
    end if
  end function intrinsic_value

  function binomial_raw(spot, strike, rate, time, sigma, dividend, payoff, steps, european) result(value)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff, steps
    logical, intent(in) :: european
    real(dp) :: value, dt, up, down, p, disc, node_spot
    real(dp), allocatable :: values(:)
    integer :: i, j
    dt = time/real(steps, dp)
    up = exp(sigma*sqrt(dt))
    down = 1.0_dp/up
    p = (exp((rate-dividend)*dt)-down)/(up-down)
    disc = exp(-rate*dt)
    allocate(values(0:steps))
    do i = 0, steps
      node_spot = spot*up**real(steps-i,dp)*down**real(i,dp)
      values(i) = intrinsic_value(node_spot, strike, payoff)
    end do
    do j = steps-1, 0, -1
      do i = 0, j
        values(i) = disc*(p*values(i)+(1.0_dp-p)*values(i+1))
        if (.not. european) then
          node_spot = spot*up**real(j-i,dp)*down**real(i,dp)
          values(i) = max(values(i), intrinsic_value(node_spot, strike, payoff))
        end if
      end do
    end do
    value = values(0)
  end function binomial_raw

  function binomial_american_price(spot, strike, rate, time, sigma, dividend, payoff, steps) result(value)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff, steps
    real(dp) :: value, amer, euro_tree, euro_exact
    amer = binomial_raw(spot,strike,rate,time,sigma,dividend,payoff,steps,.false.)
    euro_tree = binomial_raw(spot,strike,rate,time,sigma,dividend,payoff,steps,.true.)
    euro_exact = bs_european_price(spot,strike,rate,time,sigma,dividend,payoff)
    value = amer + euro_exact - euro_tree
  end function binomial_american_price

  function binomial_american_greeks(spot, strike, rate, time, sigma, dividend, payoff, steps, eps) result(res)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff
    integer, intent(in), optional :: steps
    real(dp), intent(in), optional :: eps
    type(greeks_result) :: res
    integer :: n
    real(dp) :: h, hs, hg, hv, ht, hr, hq, up, mid, down
    if (.not. ieee_is_finite(spot) .or. .not. ieee_is_finite(strike) .or. &
        .not. ieee_is_finite(rate) .or. .not. ieee_is_finite(time) .or. &
        .not. ieee_is_finite(sigma) .or. .not. ieee_is_finite(dividend) .or. &
        spot <= 0.0_dp .or. strike <= 0.0_dp .or. time <= 0.0_dp .or. sigma <= 0.0_dp) then
      res%ok = .false.
      res%message = 'spot, strike, time, and volatility must be positive finite values'
      return
    end if
    if (payoff /= payoff_call .and. payoff /= payoff_put) then
      res%ok = .false.
      res%message = 'American payoff must be call or put'
      return
    end if
    n = 1000
    if (present(steps)) n = steps
    if (n < 2) then
      res%ok = .false.
      res%message = 'steps must be at least 2'
      return
    end if
    h = 1.0e-5_dp
    if (present(eps)) h = eps
    hs = max(h, 1.0e-5_dp*spot)
    hv = min(max(h, 1.0e-6_dp), 0.25_dp*sigma)
    ht = min(max(h, 1.0e-6_dp), 0.25_dp*time)
    hr = max(h, 1.0e-6_dp)
    hq = hr
    mid = binomial_american_price(spot,strike,rate,time,sigma,dividend,payoff,n)
    res%fair_value = mid
    up = binomial_american_price(spot+hs,strike,rate,time,sigma,dividend,payoff,n)
    down = binomial_american_price(spot-hs,strike,rate,time,sigma,dividend,payoff,n)
    res%delta = (up-down)/(2.0_dp*hs)
    hg = spot/50.0_dp
    up = binomial_american_price(spot+hg,strike,rate,time,sigma,dividend,payoff,n)
    down = binomial_american_price(spot-hg,strike,rate,time,sigma,dividend,payoff,n)
    res%gamma = (up-2.0_dp*mid+down)/(hg*hg)
    up = binomial_american_price(spot,strike,rate,time,sigma+hv,dividend,payoff,n)
    down = binomial_american_price(spot,strike,rate,time,sigma-hv,dividend,payoff,n)
    res%vega = (up-down)/(2.0_dp*hv)
    up = binomial_american_price(spot,strike,rate,time+ht,sigma,dividend,payoff,n)
    down = binomial_american_price(spot,strike,rate,time-ht,sigma,dividend,payoff,n)
    res%theta = -(up-down)/(2.0_dp*ht)
    up = binomial_american_price(spot,strike,rate+hr,time,sigma,dividend,payoff,n)
    down = binomial_american_price(spot,strike,rate-hr,time,sigma,dividend,payoff,n)
    res%rho = (up-down)/(2.0_dp*hr)
    up = binomial_american_price(spot,strike,rate,time,sigma,dividend+hq,payoff,n)
    down = binomial_american_price(spot,strike,rate,time,sigma,dividend-hq,payoff,n)
    res%epsilon = (up-down)/(2.0_dp*hq)
    res%elasticity = safe_divide(spot*res%delta,res%fair_value,0.0_dp)
  end function binomial_american_greeks
end module greeks_american
