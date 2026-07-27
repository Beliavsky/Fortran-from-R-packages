! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
module greeks_implied_volatility
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use greeks_kinds, only: dp
  use greeks_types, only: implied_vol_result, greeks_result, price_callback, payoff_call, payoff_put
  use greeks_european, only: bs_european_greeks, bs_european_price
  use greeks_geometric_asian, only: bs_geometric_asian_price
  use greeks_american, only: binomial_american_price
  implicit none
  private
  public :: bs_implied_volatility, implied_volatility
  public :: geometric_asian_implied_volatility, american_implied_volatility
contains
  function bs_implied_volatility(option_price, spot, strike, rate, time, dividend, payoff, &
      start_volatility, precision, max_iter) result(res)
    real(dp), intent(in) :: option_price, spot, strike, rate, time, dividend
    integer, intent(in) :: payoff
    real(dp), intent(in), optional :: start_volatility, precision
    integer, intent(in), optional :: max_iter
    type(implied_vol_result) :: res
    real(dp) :: sigma, tol, error, denom, low, high, price_low, price_high, next
    integer :: iter, nmax
    type(greeks_result) :: g

    sigma = 0.3_dp
    if (present(start_volatility)) sigma = start_volatility
    tol = 1.0e-9_dp
    if (present(precision)) tol = precision
    nmax = 30
    if (present(max_iter)) nmax = max_iter
    if (.not. ieee_is_finite(option_price) .or. option_price < 0.0_dp .or. &
        spot <= 0.0_dp .or. strike <= 0.0_dp .or. time <= 0.0_dp .or. &
        sigma <= 0.0_dp .or. tol <= 0.0_dp .or. nmax < 1 .or. &
        (payoff /= payoff_call .and. payoff /= payoff_put)) then
      res%ok = .false.
      res%message = 'invalid implied-volatility inputs'
      return
    end if
    low = 1.0e-12_dp
    high = max(1.0_dp,2.0_dp*sigma)
    price_low = bs_european_price(spot,strike,rate,time,low,dividend,payoff)
    if (option_price <= price_low + tol) then
      res%ok = .false.
      res%message = 'option price is below the positive-volatility boundary'
      return
    end if
    price_high = bs_european_price(spot,strike,rate,time,high,dividend,payoff)
    do while (price_high < option_price .and. high < 32.0_dp)
      high = 2.0_dp*high
      price_high = bs_european_price(spot,strike,rate,time,high,dividend,payoff)
    end do
    if (price_high < option_price) then
      res%ok = .false.
      res%message = 'option price is above the search range'
      return
    end if
    sigma = min(max(sigma,low),high)
    do iter = 1, nmax
      g = bs_european_greeks(spot,strike,rate,time,sigma,dividend,payoff)
      error = g%fair_value-option_price
      if (abs(error) <= tol) then
        res%volatility = sigma
        res%price_error = error
        res%iterations = iter
        res%converged = .true.
        return
      end if
      if (error > 0.0_dp) then
        high = sigma
      else
        low = sigma
      end if
      denom = 2.0_dp*g%vega*g%vega-error*g%vomma
      if (abs(denom) > tiny(1.0_dp)) then
        next = sigma-2.0_dp*error*g%vega/denom
      else
        next = 0.5_dp*(low+high)
      end if
      if (.not. ieee_is_finite(next) .or. next <= low .or. next >= high) next = 0.5_dp*(low+high)
      sigma = next
    end do
    res%volatility = sigma
    res%price_error = bs_european_price(spot,strike,rate,time,sigma,dividend,payoff)-option_price
    res%iterations = nmax
    res%converged = abs(res%price_error) <= tol
    if (.not. res%converged) then
      res%ok = .false.
      res%message = 'maximum number of iterations reached'
    end if
  end function bs_implied_volatility

  function implied_volatility(option_price, price_fn, start_volatility, precision, max_iter) result(res)
    real(dp), intent(in) :: option_price
    procedure(price_callback) :: price_fn
    real(dp), intent(in), optional :: start_volatility, precision
    integer, intent(in), optional :: max_iter
    type(implied_vol_result) :: res
    real(dp) :: sigma, tol, low, high, p, error, h, deriv, next, plow, phigh
    integer :: iter, nmax
    sigma=0.3_dp; if (present(start_volatility)) sigma=start_volatility
    tol=1.0e-6_dp; if (present(precision)) tol=precision
    nmax=30; if (present(max_iter)) nmax=max_iter
    if (option_price < 0.0_dp .or. sigma <= 0.0_dp .or. tol <= 0.0_dp .or. nmax < 1) then
      res%ok=.false.; res%message='invalid implied-volatility inputs'; return
    end if
    low=1.0e-8_dp; high=max(1.0_dp,2.0_dp*sigma)
    plow=price_fn(low); phigh=price_fn(high)
    do while (phigh < option_price .and. high < 32.0_dp)
      high=2.0_dp*high; phigh=price_fn(high)
    end do
    if (.not. ieee_is_finite(plow) .or. .not. ieee_is_finite(phigh) .or. &
        option_price < plow-tol .or. option_price > phigh+tol) then
      res%ok=.false.; res%message='target price is outside the volatility bracket'; return
    end if
    sigma=min(max(sigma,low),high)
    do iter=1,nmax
      p=price_fn(sigma); error=p-option_price
      if (abs(error) <= tol) then
        res%volatility=sigma; res%price_error=error; res%iterations=iter; res%converged=.true.; return
      end if
      if (error > 0.0_dp) then; high=sigma; else; low=sigma; end if
      h=max(1.0e-5_dp,1.0e-4_dp*sigma)
      deriv=(price_fn(sigma+h)-price_fn(max(low,sigma-h)))/(sigma+h-max(low,sigma-h))
      if (deriv > tiny(1.0_dp)) then; next=sigma-error/deriv; else; next=0.5_dp*(low+high); end if
      if (.not. ieee_is_finite(next) .or. next <= low .or. next >= high) next=0.5_dp*(low+high)
      sigma=next
    end do
    res%volatility=sigma; res%price_error=price_fn(sigma)-option_price; res%iterations=nmax
    res%converged=abs(res%price_error)<=tol
    if (.not. res%converged) then; res%ok=.false.; res%message='maximum number of iterations reached'; end if
  end function implied_volatility

  function builtin_price(style, volatility, spot, strike, rate, time, dividend, payoff, steps) result(value)
    integer, intent(in) :: style, payoff, steps
    real(dp), intent(in) :: volatility, spot, strike, rate, time, dividend
    real(dp) :: value
    if (style == 1) then
      value=bs_geometric_asian_price(spot,strike,rate,time,volatility,dividend,payoff)
    else
      value=binomial_american_price(spot,strike,rate,time,volatility,dividend,payoff,steps)
    end if
  end function builtin_price

  function builtin_implied_volatility(style, option_price, spot, strike, rate, time, dividend, payoff, &
      steps, start, tol, nmax) result(res)
    integer, intent(in) :: style, payoff, steps, nmax
    real(dp), intent(in) :: option_price, spot, strike, rate, time, dividend, start, tol
    type(implied_vol_result) :: res
    real(dp) :: sigma, low, high, plow, phigh, p, error, h, deriv, next, lower_eval
    integer :: iter
    if (option_price < 0.0_dp .or. start <= 0.0_dp .or. tol <= 0.0_dp .or. nmax < 1) then
      res%ok=.false.; res%message='invalid implied-volatility inputs'; return
    end if
    low=3.0e-3_dp; high=max(1.0_dp,2.0_dp*start)
    plow=builtin_price(style,low,spot,strike,rate,time,dividend,payoff,steps)
    phigh=builtin_price(style,high,spot,strike,rate,time,dividend,payoff,steps)
    do while (phigh < option_price .and. high < 32.0_dp)
      high=2.0_dp*high
      phigh=builtin_price(style,high,spot,strike,rate,time,dividend,payoff,steps)
    end do
    if (.not. ieee_is_finite(plow) .or. .not. ieee_is_finite(phigh) .or. &
        option_price < plow-tol .or. option_price > phigh+tol) then
      res%ok=.false.; res%message='target price is outside the volatility bracket'; return
    end if
    sigma=min(max(start,low),high)
    do iter=1,nmax
      p=builtin_price(style,sigma,spot,strike,rate,time,dividend,payoff,steps)
      error=p-option_price
      if (abs(error) <= tol) then
        res%volatility=sigma; res%price_error=error; res%iterations=iter; res%converged=.true.; return
      end if
      if (error > 0.0_dp) then; high=sigma; else; low=sigma; end if
      h=max(1.0e-5_dp,1.0e-4_dp*sigma)
      lower_eval=max(low,sigma-h)
      deriv=(builtin_price(style,sigma+h,spot,strike,rate,time,dividend,payoff,steps)- &
        builtin_price(style,lower_eval,spot,strike,rate,time,dividend,payoff,steps))/(sigma+h-lower_eval)
      if (deriv > tiny(1.0_dp)) then; next=sigma-error/deriv; else; next=0.5_dp*(low+high); end if
      if (.not. ieee_is_finite(next) .or. next <= low .or. next >= high) next=0.5_dp*(low+high)
      sigma=next
    end do
    res%volatility=sigma
    res%price_error=builtin_price(style,sigma,spot,strike,rate,time,dividend,payoff,steps)-option_price
    res%iterations=nmax; res%converged=abs(res%price_error)<=tol
    if (.not. res%converged) then; res%ok=.false.; res%message='maximum number of iterations reached'; end if
  end function builtin_implied_volatility

  function geometric_asian_implied_volatility(option_price, spot, strike, rate, time, dividend, payoff, &
      start_volatility, precision, max_iter) result(res)
    real(dp), intent(in) :: option_price, spot, strike, rate, time, dividend
    integer, intent(in) :: payoff
    real(dp), intent(in), optional :: start_volatility, precision
    integer, intent(in), optional :: max_iter
    type(implied_vol_result) :: res
    real(dp) :: start, tol
    integer :: nmax
    start=0.3_dp; if (present(start_volatility)) start=start_volatility
    tol=1.0e-6_dp; if (present(precision)) tol=precision
    nmax=30; if (present(max_iter)) nmax=max_iter
    res=builtin_implied_volatility(1,option_price,spot,strike,rate,time,dividend,payoff,1,start,tol,nmax)
  end function geometric_asian_implied_volatility

  function american_implied_volatility(option_price, spot, strike, rate, time, dividend, payoff, steps, &
      start_volatility, precision, max_iter) result(res)
    real(dp), intent(in) :: option_price, spot, strike, rate, time, dividend
    integer, intent(in) :: payoff
    integer, intent(in), optional :: steps, max_iter
    real(dp), intent(in), optional :: start_volatility, precision
    type(implied_vol_result) :: res
    integer :: n, nmax
    real(dp) :: start, tol
    n=500; if (present(steps)) n=steps
    start=0.3_dp; if (present(start_volatility)) start=start_volatility
    tol=1.0e-6_dp; if (present(precision)) tol=precision
    nmax=30; if (present(max_iter)) nmax=max_iter
    res=builtin_implied_volatility(2,option_price,spot,strike,rate,time,dividend,payoff,n,start,tol,nmax)
  end function american_implied_volatility
end module greeks_implied_volatility
