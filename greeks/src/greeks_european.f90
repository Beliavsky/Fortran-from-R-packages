! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
module greeks_european
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use greeks_kinds, only: dp
  use greeks_types, only: greeks_result, payoff_call, payoff_put, &
    payoff_cash_call, payoff_cash_put, payoff_asset_call, payoff_asset_put
  use greeks_math, only: normal_pdf, normal_cdf, safe_divide
  implicit none
  private
  public :: bs_european_greeks, bs_european_price
contains
  pure function bs_european_price(spot, strike, rate, time, sigma, dividend, payoff) result(value)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff
    real(dp) :: value, d1, d2, ds, dk
    d1 = (log(spot/strike) + (rate-dividend+0.5_dp*sigma*sigma)*time)/(sigma*sqrt(time))
    d2 = d1 - sigma*sqrt(time)
    ds = spot*exp(-dividend*time)
    dk = strike*exp(-rate*time)
    select case (payoff)
    case (payoff_call)
      value = ds*normal_cdf(d1) - dk*normal_cdf(d2)
    case (payoff_put)
      value = dk*normal_cdf(-d2) - ds*normal_cdf(-d1)
    case (payoff_cash_call)
      value = exp(-rate*time)*normal_cdf(d2)
    case (payoff_cash_put)
      value = exp(-rate*time)*normal_cdf(-d2)
    case (payoff_asset_call)
      value = ds*normal_cdf(d1)
    case (payoff_asset_put)
      value = ds*normal_cdf(-d1)
    case default
      value = 0.0_dp
    end select
  end function bs_european_price

  pure function first_derivative(prefactor, cdfx, pdfx, log_i, x_i) result(value)
    real(dp), intent(in) :: prefactor, cdfx, pdfx, log_i, x_i
    real(dp) :: value
    value = prefactor*(log_i*cdfx + pdfx*x_i)
  end function first_derivative

  pure function second_derivative(prefactor, cdfx, pdfx, x, log_i, log_j, &
      log_ij, x_i, x_j, x_ij) result(value)
    real(dp), intent(in) :: prefactor, cdfx, pdfx, x
    real(dp), intent(in) :: log_i, log_j, log_ij, x_i, x_j, x_ij
    real(dp) :: value, a, b, c
    a = log_ij + log_i*log_j
    b = log_i*x_j + log_j*x_i + x_ij
    c = x_i*x_j
    value = prefactor*(a*cdfx + b*pdfx - c*x*pdfx)
  end function second_derivative

  pure function third_derivative(prefactor, cdfx, pdfx, x, &
      log_i, log_j, log_k, log_ij, log_ik, log_jk, log_ijk, &
      x_i, x_j, x_k, x_ij, x_ik, x_jk, x_ijk) result(value)
    real(dp), intent(in) :: prefactor, cdfx, pdfx, x
    real(dp), intent(in) :: log_i, log_j, log_k, log_ij, log_ik, log_jk, log_ijk
    real(dp), intent(in) :: x_i, x_j, x_k, x_ij, x_ik, x_jk, x_ijk
    real(dp) :: value, a, b, c, ak, bk, ck
    a = log_ij + log_i*log_j
    b = log_i*x_j + log_j*x_i + x_ij
    c = x_i*x_j
    ak = log_ijk + log_ik*log_j + log_i*log_jk
    bk = log_ik*x_j + log_i*x_jk + log_jk*x_i + log_j*x_ik + x_ijk
    ck = x_ik*x_j + x_i*x_jk
    value = prefactor*((log_k*a+ak)*cdfx + (log_k*b+a*x_k+bk)*pdfx - &
      (log_k*c+b*x_k+ck)*x*pdfx + c*x_k*(x*x-1.0_dp)*pdfx)
  end function third_derivative

  function bs_european_greeks(spot, strike, rate, time, sigma, dividend, payoff) result(res)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    integer, intent(in) :: payoff
    type(greeks_result) :: res
    real(dp) :: sqrt_t, d1, d2, eq, er, nd1, pd1, pd2, pm1, pm2, value
    logical :: asset, is_call
    real(dp) :: signx, shift, base_d, x, cdfx, pdfx, pref
    real(dp) :: lp_s, lp_v, lp_t, lp_r, lp_q, lp_ss, lp_sss
    real(dp) :: d0, ds, dv, dt, dr, dq, dss, dsv, dst, dvv, dvt, dvr
    real(dp) :: dsss, dsst, dssv, dvvv
    real(dp) :: xs, xv, xt, xr, xq, xss, xsv, xst, xvv, xvt, xvr
    real(dp) :: xsss, xsst, xssv, xvvv

    if (.not. ieee_is_finite(spot) .or. .not. ieee_is_finite(strike) .or. &
        .not. ieee_is_finite(rate) .or. .not. ieee_is_finite(time) .or. &
        .not. ieee_is_finite(sigma) .or. .not. ieee_is_finite(dividend) .or. &
        spot <= 0.0_dp .or. strike <= 0.0_dp .or. time <= 0.0_dp .or. sigma <= 0.0_dp) then
      res%ok = .false.
      res%message = 'spot, strike, time, and volatility must be positive finite values'
      return
    end if
    if (payoff < payoff_call .or. payoff > payoff_asset_put) then
      res%ok = .false.
      res%message = 'unsupported payoff code'
      return
    end if

    sqrt_t = sqrt(time)
    d1 = (log(spot/strike) + (rate-dividend+0.5_dp*sigma*sigma)*time)/(sigma*sqrt_t)
    d2 = d1 - sigma*sqrt_t
    eq = exp(-dividend*time)
    er = exp(-rate*time)
    nd1 = normal_pdf(d1)
    pd1 = normal_cdf(d1)
    pd2 = normal_cdf(d2)

    if (payoff == payoff_call .or. payoff == payoff_put) then
      if (payoff == payoff_call) then
        value = spot*eq*pd1 - strike*er*pd2
        res%fair_value = value
        res%delta = eq*pd1
        res%theta = -spot*nd1*sigma*eq/(2.0_dp*sqrt_t) + dividend*spot*pd1*eq - rate*strike*er*pd2
        res%rho = strike*time*er*pd2
        res%epsilon = -spot*time*eq*pd1
        res%charm = dividend*eq*pd1 - eq*nd1*(2.0_dp*(rate-dividend)*time-d2*sigma*sqrt_t) / &
          (2.0_dp*time*sigma*sqrt_t)
      else
        pm1 = normal_cdf(-d1)
        pm2 = normal_cdf(-d2)
        value = strike*er*pm2 - spot*eq*pm1
        res%fair_value = value
        res%delta = -eq*pm1
        res%theta = -spot*nd1*sigma*eq/(2.0_dp*sqrt_t) - dividend*spot*pm1*eq + rate*strike*er*pm2
        res%rho = -strike*time*er*pm2
        res%epsilon = spot*time*eq*pm1
        res%charm = -dividend*eq*pm1 - eq*nd1*(2.0_dp*(rate-dividend)*time-d2*sigma*sqrt_t) / &
          (2.0_dp*time*sigma*sqrt_t)
      end if
      res%elasticity = safe_divide(spot*res%delta, value, 0.0_dp)
      res%vega = spot*eq*nd1*sqrt_t
      res%gamma = nd1*eq/(spot*sigma*sqrt_t)
      res%vanna = -eq*nd1*d2/sigma
      res%vomma = res%vega*d1*d2/sigma
      res%veta = -spot*eq*nd1*sqrt_t*(dividend + (rate-dividend)*d1/(sigma*sqrt_t) - &
        (1.0_dp+d1*d2)/(2.0_dp*time))
      res%vera = strike*time*er*normal_pdf(d2)*(-sqrt_t-d2/sigma)
      res%speed = -eq*nd1/(spot*spot*sigma*sqrt_t)*(d1/(sigma*sqrt_t)+1.0_dp)
      res%zomma = eq*nd1*(d1*d2-1.0_dp)/(spot*sigma*sigma*sqrt_t)
      res%color = -eq*nd1/(2.0_dp*spot*time**1.5_dp*sigma) * &
        (2.0_dp*dividend*time+1.0_dp + &
         (2.0_dp*(rate-dividend)*time-d2*sigma*sqrt_t)/(sigma*sqrt_t)*d1)
      res%ultima = -res%vega/(sigma*sigma)*(d1*d2*(1.0_dp-d1*d2)+d1*d1+d2*d2)
      return
    end if

    asset = payoff == payoff_asset_call .or. payoff == payoff_asset_put
    is_call = payoff == payoff_cash_call .or. payoff == payoff_asset_call
    if (is_call) then
      signx = 1.0_dp
    else
      signx = -1.0_dp
    end if
    if (asset) then
      shift = 0.5_dp
      base_d = d1
      pref = spot*eq
      lp_s = 1.0_dp/spot
      lp_t = -dividend
      lp_r = 0.0_dp
      lp_q = -time
      lp_ss = -1.0_dp/(spot*spot)
      lp_sss = 2.0_dp/(spot**3)
    else
      shift = -0.5_dp
      base_d = d2
      pref = er
      lp_s = 0.0_dp
      lp_t = -rate
      lp_r = -time
      lp_q = 0.0_dp
      lp_ss = 0.0_dp
      lp_sss = 0.0_dp
    end if
    lp_v = 0.0_dp
    x = signx*base_d
    cdfx = normal_cdf(x)
    pdfx = normal_pdf(x)
    d0 = base_d - shift*sigma*sqrt_t
    ds = 1.0_dp/(spot*sigma*sqrt_t)
    dv = -base_d/sigma + 2.0_dp*shift*sqrt_t
    dt = (rate-dividend+shift*sigma*sigma)/(sigma*sqrt_t) - base_d/(2.0_dp*time)
    dr = sqrt_t/sigma
    dq = -sqrt_t/sigma
    dss = -1.0_dp/(spot*spot*sigma*sqrt_t)
    dsv = -1.0_dp/(spot*sigma*sigma*sqrt_t)
    dst = -1.0_dp/(2.0_dp*spot*sigma*time*sqrt_t)
    dvv = 2.0_dp*d0/(sigma*sigma)
    dvt = -dt/sigma + shift/sqrt_t
    dvr = -sqrt_t/(sigma*sigma)
    dsss = 2.0_dp/(spot**3*sigma*sqrt_t)
    dsst = 1.0_dp/(2.0_dp*spot*spot*sigma*time*sqrt_t)
    dssv = 1.0_dp/(spot*spot*sigma*sigma*sqrt_t)
    dvvv = -6.0_dp*d0/(sigma**3)
    xs=signx*ds; xv=signx*dv; xt=signx*dt; xr=signx*dr; xq=signx*dq
    xss=signx*dss; xsv=signx*dsv; xst=signx*dst; xvv=signx*dvv
    xvt=signx*dvt; xvr=signx*dvr
    xsss=signx*dsss; xsst=signx*dsst; xssv=signx*dssv; xvvv=signx*dvvv

    res%fair_value = pref*cdfx
    res%delta = first_derivative(pref,cdfx,pdfx,lp_s,xs)
    res%vega = first_derivative(pref,cdfx,pdfx,lp_v,xv)
    res%theta = -first_derivative(pref,cdfx,pdfx,lp_t,xt)
    res%rho = first_derivative(pref,cdfx,pdfx,lp_r,xr)
    res%epsilon = first_derivative(pref,cdfx,pdfx,lp_q,xq)
    res%elasticity = safe_divide(res%delta*spot,res%fair_value,0.0_dp)
    res%gamma = second_derivative(pref,cdfx,pdfx,x,lp_s,lp_s,lp_ss,xs,xs,xss)
    res%vanna = second_derivative(pref,cdfx,pdfx,x,lp_s,lp_v,0.0_dp,xs,xv,xsv)
    res%charm = -second_derivative(pref,cdfx,pdfx,x,lp_s,lp_t,0.0_dp,xs,xt,xst)
    res%vomma = second_derivative(pref,cdfx,pdfx,x,lp_v,lp_v,0.0_dp,xv,xv,xvv)
    res%veta = second_derivative(pref,cdfx,pdfx,x,lp_v,lp_t,0.0_dp,xv,xt,xvt)
    res%vera = second_derivative(pref,cdfx,pdfx,x,lp_r,lp_v,0.0_dp,xr,xv,xvr)
    res%speed = third_derivative(pref,cdfx,pdfx,x,lp_s,lp_s,lp_s,lp_ss,lp_ss,lp_ss,lp_sss, &
      xs,xs,xs,xss,xss,xss,xsss)
    res%zomma = third_derivative(pref,cdfx,pdfx,x,lp_s,lp_s,lp_v,lp_ss,0.0_dp,0.0_dp,0.0_dp, &
      xs,xs,xv,xss,xsv,xsv,xssv)
    res%color = third_derivative(pref,cdfx,pdfx,x,lp_s,lp_s,lp_t,lp_ss,0.0_dp,0.0_dp,0.0_dp, &
      xs,xs,xt,xss,xst,xst,xsst)
    res%ultima = third_derivative(pref,cdfx,pdfx,x,lp_v,lp_v,lp_v,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
      xv,xv,xv,xvv,xvv,xvv,xvvv)
  end function bs_european_greeks
end module greeks_european
