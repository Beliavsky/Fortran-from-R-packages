! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_black_scholes
  use qrmtools_kinds, only : dp
  use qrmtools_types, only : greeks_result
  use qrmtools_stats, only : normal_pdf, normal_cdf
  implicit none
  private
  public :: black_scholes, black_scholes_greeks
contains
  pure real(dp) function black_scholes(t,spot,rate,volatility,strike,maturity,is_put) result(price)
    real(dp), intent(in) :: t,spot,rate,volatility,strike,maturity
    logical, intent(in), optional :: is_put
    real(dp) :: tau,d1,d2
    logical :: put
    put=.false.; if(present(is_put))put=is_put
    tau=maturity-t
    if(tau<=0.0_dp) then
      if(put) then; price=max(strike-spot,0.0_dp); else; price=max(spot-strike,0.0_dp); end if
      return
    end if
    if(volatility<=0.0_dp) then
      if(put) then; price=max(strike*exp(-rate*tau)-spot,0.0_dp)
      else; price=max(spot-strike*exp(-rate*tau),0.0_dp); end if
      return
    end if
    d1=(log(spot/strike)+(rate+0.5_dp*volatility**2)*tau)/(volatility*sqrt(tau))
    d2=d1-volatility*sqrt(tau)
    if(put) then
      price=-spot*normal_cdf(-d1)+strike*exp(-rate*tau)*normal_cdf(-d2)
    else
      price=spot*normal_cdf(d1)-strike*exp(-rate*tau)*normal_cdf(d2)
    end if
  end function black_scholes

  pure function black_scholes_greeks(t,spot,rate,volatility,strike,maturity,is_put) result(output)
    real(dp), intent(in) :: t,spot,rate,volatility,strike,maturity
    logical, intent(in), optional :: is_put
    type(greeks_result) :: output
    real(dp) :: tau,d1,d2
    logical :: put
    put=.false.; if(present(is_put))put=is_put
    tau=maturity-t
    if(tau<=0.0_dp .or. spot<=0.0_dp .or. strike<=0.0_dp .or. volatility<=0.0_dp) then
      output%message='Positive time, spot, strike, and volatility are required.'; return
    end if
    d1=(log(spot/strike)+(rate+0.5_dp*volatility**2)*tau)/(volatility*sqrt(tau))
    d2=d1-volatility*sqrt(tau)
    output%price=black_scholes(t,spot,rate,volatility,strike,maturity,put)
    output%vega=spot*normal_pdf(d1)*sqrt(tau)
    output%gamma=normal_pdf(d1)/(spot*volatility*sqrt(tau))
    output%vanna=-normal_pdf(d1)*d2/volatility
    output%vomma=output%vega*d1*d2/volatility
    if(put) then
      output%delta=-normal_cdf(-d1)
      output%theta=-(spot*normal_pdf(d1)*volatility)/(2.0_dp*sqrt(tau))+&
        rate*strike*exp(-rate*tau)*normal_cdf(-d2)
      output%rho=-strike*tau*exp(-rate*tau)*normal_cdf(-d2)
    else
      output%delta=normal_cdf(d1)
      output%theta=-(spot*normal_pdf(d1)*volatility)/(2.0_dp*sqrt(tau))-&
        rate*strike*exp(-rate*tau)*normal_cdf(d2)
      output%rho=strike*tau*exp(-rate*tau)*normal_cdf(d2)
    end if
    output%ok=.true.
  end function black_scholes_greeks
end module qrmtools_black_scholes
