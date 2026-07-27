! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
module greeks_geometric_asian
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use greeks_kinds, only: dp
  use greeks_types, only: greeks_result, payoff_call, payoff_put
  use greeks_math, only: normal_pdf, normal_cdf, safe_divide
  implicit none
  private
  public :: bs_geometric_asian_greeks, bs_geometric_asian_price
contains
  pure function bs_geometric_asian_price(spot, strike, rate, time, sigma, dividend, payoff) result(value)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff
    real(dp) :: value, d, shift, pref, dk
    d = (log(spot/strike) + 0.5_dp*time*(rate-dividend-0.5_dp*sigma*sigma)) / &
      (sigma*sqrt(time/3.0_dp))
    shift = sigma*sqrt(time/3.0_dp)
    pref = exp(-0.5_dp*time*(rate+dividend+sigma*sigma/6.0_dp))
    dk = strike*exp(-rate*time)
    if (payoff == payoff_call) then
      value = spot*pref*normal_cdf(d+shift) - dk*normal_cdf(d)
    else if (payoff == payoff_put) then
      value = -spot*pref*normal_cdf(-d-shift) + dk*normal_cdf(-d)
    else
      value = 0.0_dp
    end if
  end function bs_geometric_asian_price

  function bs_geometric_asian_greeks(spot, strike, rate, time, sigma, dividend, payoff) result(res)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff
    type(greeks_result) :: res
    real(dp) :: rmq, lm, st, st3, s3, s3t, d, dd, dv, dt, dg
    real(dp) :: shift, ad, disc, pref, dk, shift_v, dvv, adv, advv, pvr, pvrr
    real(dp) :: pa, ma

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
      res%message = 'geometric Asian payoff must be call or put'
      return
    end if

    rmq = rate-dividend
    lm = log(spot/strike)
    st = sqrt(time)
    st3 = sqrt(time/3.0_dp)
    s3 = sqrt(3.0_dp)
    s3t = sqrt(3.0_dp*time)
    d = (lm+0.5_dp*time*(rmq-0.5_dp*sigma*sigma))/(sigma*st3)
    dd = s3/(spot*sigma*st)
    dv = -s3t*(0.25_dp + lm/(sigma*sigma*time) + rmq/(2.0_dp*sigma*sigma))
    dt = -s3*lm/(2.0_dp*sigma*time**1.5_dp) + s3*(rmq-0.5_dp*sigma*sigma)/(4.0_dp*sigma*st)
    dg = -s3/(spot*spot*sigma*st)
    shift = sigma*st3
    ad = d+shift
    disc = exp(-rate*time)
    pref = exp(-0.5_dp*time*(rate+dividend+sigma*sigma/6.0_dp))
    dk = disc*strike
    shift_v = st3
    dvv = 2.0_dp*s3t*(lm/time+rmq/2.0_dp)/(sigma**3)
    adv = dv+shift_v
    advv = dvv
    pvr = -sigma*time/6.0_dp
    pvrr = -time/6.0_dp

    if (payoff == payoff_call) then
      res%fair_value = spot*pref*normal_cdf(ad)-dk*normal_cdf(d)
      res%delta = pref*(normal_cdf(ad)+spot*normal_pdf(ad)*dd)-dk*normal_pdf(d)*dd
      res%vega = spot*pref*(pvr*normal_cdf(ad)+normal_pdf(ad)*adv)-dk*normal_pdf(d)*dv
      res%rho = spot*pref*(-0.5_dp*time*normal_cdf(ad)+normal_pdf(ad)*s3t/(2.0_dp*sigma)) + &
        strike*disc*(time*normal_cdf(d)-normal_pdf(d)*s3t/(2.0_dp*sigma))
      res%theta = -spot*pref*(-0.5_dp*(rate+dividend+sigma*sigma/6.0_dp)*normal_cdf(ad) + &
        normal_pdf(ad)*(dt+sigma/(2.0_dp*sqrt(3.0_dp*time)))) + &
        strike*disc*(-rate*normal_cdf(d)+normal_pdf(d)*dt)
      res%gamma = pref*normal_pdf(ad)*(2.0_dp*dd-spot*ad*dd*dd+spot*dg) + &
        dk*normal_pdf(d)*(d*dd*dd-dg)
      res%vomma = spot*pref*((pvr*pvr+pvrr)*normal_cdf(ad) + &
        (2.0_dp*pvr*adv-ad*adv*adv+advv)*normal_pdf(ad)) + &
        dk*normal_pdf(d)*(d*dv*dv-dvv)
    else
      pa = -ad
      ma = -d
      res%fair_value = -spot*pref*normal_cdf(pa)+dk*normal_cdf(ma)
      res%delta = -pref*(normal_cdf(pa)-spot*normal_pdf(pa)*dd)-dk*normal_pdf(ma)*dd
      res%vega = -spot*pref*(pvr*normal_cdf(pa)-normal_pdf(pa)*adv)-dk*normal_pdf(ma)*dv
      res%rho = -spot*pref*(-0.5_dp*time*normal_cdf(pa)-normal_pdf(pa)*s3t/(2.0_dp*sigma)) - &
        strike*disc*(time*normal_cdf(ma)+normal_pdf(ma)*s3t/(2.0_dp*sigma))
      res%theta = spot*pref*(-0.5_dp*(rate+dividend+sigma*sigma/6.0_dp)*normal_cdf(pa) - &
        normal_pdf(pa)*(dt+sigma/(2.0_dp*sqrt(3.0_dp*time)))) + &
        strike*disc*(rate*normal_cdf(ma)+normal_pdf(ma)*dt)
      res%gamma = pref*normal_pdf(pa)*(2.0_dp*dd-spot*pa*dd*dd+spot*dg) - &
        dk*normal_pdf(ma)*(d*dd*dd-dg)
      res%vomma = -spot*pref*((pvr*pvr+pvrr)*normal_cdf(pa) + &
        (-2.0_dp*pvr*adv+ad*adv*adv-advv)*normal_pdf(ad)) + &
        dk*normal_pdf(d)*(d*dv*dv-dvv)
    end if
    res%elasticity = safe_divide(spot*res%delta, res%fair_value, 0.0_dp)
  end function bs_geometric_asian_greeks
end module greeks_geometric_asian
