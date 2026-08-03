! SPDX-License-Identifier: MIT
module greeks_black_scholes
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use greeks_kinds, only: dp
  use greeks_math, only: normal_cdf, normal_pdf
  use greeks_types, only: greek_result, initialize_result, set_error
  use greeks_types, only: greeks_invalid_argument, greeks_unknown_payoff
  use greeks_types, only: greeks_no_convergence, greeks_numerical_error
  implicit none
  private
  public :: bs_european_greeks, bs_geometric_asian_greeks
  public :: bs_implied_volatility, bs_european_price

contains

  subroutine validate_inputs(spot, strike, time, sigma, result, valid)
    real(dp), intent(in) :: spot, strike, time, sigma
    type(greek_result), intent(inout) :: result
    logical, intent(out) :: valid
    valid = .false.
    if (spot <= 0.0_dp .or. strike <= 0.0_dp .or. time <= 0.0_dp .or. &
        sigma <= 0.0_dp) then
      call set_error(result, greeks_invalid_argument, &
        'spot, strike, time, and volatility must be positive')
      return
    end if
    valid = .true.
  end subroutine validate_inputs

  pure real(dp) function bs_european_price(spot, strike, rate, time, sigma, &
      dividend, payoff) result(value)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    character(len=*), intent(in) :: payoff
    real(dp) :: d1, d2, sqrt_time, ds, dk

    sqrt_time = sqrt(time)
    d1 = (log(spot/strike) + (rate - dividend + 0.5_dp*sigma**2)*time) / &
      (sigma*sqrt_time)
    d2 = d1 - sigma*sqrt_time
    ds = spot*exp(-dividend*time)
    dk = strike*exp(-rate*time)
    select case (trim(payoff))
    case ('call')
      value = ds*normal_cdf(d1) - dk*normal_cdf(d2)
    case ('put')
      value = dk*normal_cdf(-d2) - ds*normal_cdf(-d1)
    case ('cash_or_nothing_call', 'digital_call')
      value = exp(-rate*time)*normal_cdf(d2)
    case ('cash_or_nothing_put', 'digital_put')
      value = exp(-rate*time)*normal_cdf(-d2)
    case ('asset_or_nothing_call')
      value = ds*normal_cdf(d1)
    case ('asset_or_nothing_put')
      value = ds*normal_cdf(-d1)
    case default
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    end select
  end function bs_european_price

  subroutine bs_european_greeks(spot, strike, rate, time, sigma, dividend, &
      payoff, requested, result)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    character(len=*), intent(in) :: payoff
    character(len=*), intent(in) :: requested(:)
    type(greek_result), intent(out) :: result
    real(dp) :: d1, d2, sqrt_time, eqt, ert, nd1, nd2, pd1
    real(dp) :: fair_value, value
    integer :: i
    logical :: valid

    call initialize_result(result, requested)
    call validate_inputs(spot, strike, time, sigma, result, valid)
    if (.not. valid) return
    if (.not. valid_european_payoff(payoff)) then
      call set_error(result, greeks_unknown_payoff, 'unsupported European payoff')
      return
    end if

    sqrt_time = sqrt(time)
    d1 = (log(spot/strike) + (rate - dividend + 0.5_dp*sigma**2)*time) / &
      (sigma*sqrt_time)
    d2 = d1 - sigma*sqrt_time
    eqt = exp(-dividend*time)
    ert = exp(-rate*time)
    nd1 = normal_cdf(d1)
    nd2 = normal_cdf(d2)
    pd1 = normal_pdf(d1)
    fair_value = bs_european_price(spot, strike, rate, time, sigma, dividend, payoff)

    do i = 1, size(requested)
      if (trim(payoff) == 'call' .or. trim(payoff) == 'put') then
        value = vanilla_greek(trim(requested(i)), trim(payoff), spot, strike, &
          rate, time, sigma, dividend, d1, d2, eqt, ert, nd1, nd2, pd1, fair_value)
      else
        value = binary_greek(trim(requested(i)), trim(payoff), spot, rate, time, &
          sigma, dividend, d1, d2, fair_value)
      end if
      if (.not. ieee_is_finite_local(value)) then
        call set_error(result, greeks_numerical_error, &
          'unknown Greek name or nonfinite formula result')
        return
      end if
      result%values(i) = value
    end do
  end subroutine bs_european_greeks

  pure logical function valid_european_payoff(payoff) result(valid)
    character(len=*), intent(in) :: payoff
    select case (trim(payoff))
    case ('call', 'put', 'cash_or_nothing_call', 'cash_or_nothing_put', &
          'asset_or_nothing_call', 'asset_or_nothing_put', &
          'digital_call', 'digital_put')
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_european_payoff

  pure logical function ieee_is_finite_local(x) result(value)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(dp), intent(in) :: x
    value = ieee_is_finite(x)
  end function ieee_is_finite_local

  pure real(dp) function vanilla_greek(name, payoff, spot, strike, rate, time, &
      sigma, dividend, d1, d2, eqt, ert, nd1, nd2, pd1, fair_value) result(value)
    character(len=*), intent(in) :: name, payoff
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    real(dp), intent(in) :: d1, d2, eqt, ert, nd1, nd2, pd1, fair_value
    real(dp) :: sqrt_time, nmd1, nmd2

    sqrt_time = sqrt(time)
    nmd1 = normal_cdf(-d1)
    nmd2 = normal_cdf(-d2)
    value = ieee_value(0.0_dp, ieee_quiet_nan)
    select case (name)
    case ('fair_value')
      value = fair_value
    case ('delta')
      if (payoff == 'call') then
        value = eqt*nd1
      else
        value = -eqt*nmd1
      end if
    case ('vega')
      value = spot*eqt*pd1*sqrt_time
    case ('theta')
      if (payoff == 'call') then
        value = -spot*pd1*sigma*eqt/(2.0_dp*sqrt_time) + &
          dividend*spot*nd1*eqt - rate*strike*ert*nd2
      else
        value = -spot*pd1*sigma*eqt/(2.0_dp*sqrt_time) - &
          dividend*spot*nmd1*eqt + rate*strike*ert*nmd2
      end if
    case ('rho')
      if (payoff == 'call') then
        value = strike*time*ert*nd2
      else
        value = -strike*time*ert*nmd2
      end if
    case ('epsilon')
      if (payoff == 'call') then
        value = -spot*time*eqt*nd1
      else
        value = spot*time*eqt*nmd1
      end if
    case ('lambda')
      if (payoff == 'call') then
        value = spot*eqt*nd1/fair_value
      else
        value = -spot*eqt*nmd1/fair_value
      end if
    case ('gamma')
      value = pd1*eqt/(spot*sigma*sqrt_time)
    case ('vanna')
      value = -eqt*pd1*d2/sigma
    case ('charm')
      if (payoff == 'call') then
        value = dividend*eqt*nd1 - eqt*pd1*(2.0_dp*(rate - dividend)*time - &
          d2*sigma*sqrt_time)/(2.0_dp*time*sigma*sqrt_time)
      else
        value = -dividend*eqt*nmd1 - eqt*pd1*(2.0_dp*(rate - dividend)*time - &
          d2*sigma*sqrt_time)/(2.0_dp*time*sigma*sqrt_time)
      end if
    case ('vomma')
      value = spot*eqt*pd1*sqrt_time*d1*d2/sigma
    case ('veta')
      value = -spot*eqt*pd1*sqrt_time*(dividend + &
        (rate - dividend)*d1/(sigma*sqrt_time) - (1.0_dp + d1*d2)/(2.0_dp*time))
    case ('vera')
      value = strike*time*ert*normal_pdf(d2)*(-sqrt_time - d2/sigma)
    case ('speed')
      value = -eqt*pd1/(spot**2*sigma*sqrt_time) * &
        (d1/(sigma*sqrt_time) + 1.0_dp)
    case ('zomma')
      value = eqt*pd1*(d1*d2 - 1.0_dp)/(spot*sigma**2*sqrt_time)
    case ('color')
      value = -eqt*pd1/(2.0_dp*spot*time**1.5_dp*sigma) * &
        (2.0_dp*dividend*time + 1.0_dp + &
        ((2.0_dp*(rate - dividend)*time - d2*sigma*sqrt_time) / &
        (sigma*sqrt_time))*d1)
    case ('ultima')
      value = -spot*eqt*pd1*sqrt_time/sigma**2 * &
        (d1*d2*(1.0_dp - d1*d2) + d1**2 + d2**2)
    end select
  end function vanilla_greek

  pure real(dp) function first_derivative(prefactor, cdf_value, pdf_value, &
      log_i, x_i) result(value)
    real(dp), intent(in) :: prefactor, cdf_value, pdf_value, log_i, x_i
    value = prefactor*(log_i*cdf_value + pdf_value*x_i)
  end function first_derivative

  pure real(dp) function second_derivative(prefactor, cdf_value, pdf_value, x, &
      log_i, log_j, log_ij, x_i, x_j, x_ij) result(value)
    real(dp), intent(in) :: prefactor, cdf_value, pdf_value, x
    real(dp), intent(in) :: log_i, log_j, log_ij, x_i, x_j, x_ij
    value = prefactor*((log_ij + log_i*log_j)*cdf_value + &
      (log_i*x_j + log_j*x_i + x_ij)*pdf_value - x_i*x_j*x*pdf_value)
  end function second_derivative

  pure real(dp) function third_derivative(prefactor, cdf_value, pdf_value, x, &
      li, lj, lk, lij, lik, ljk, lijk, xi, xj, xk, xij, xik, xjk, xijk) &
      result(value)
    real(dp), intent(in) :: prefactor, cdf_value, pdf_value, x
    real(dp), intent(in) :: li, lj, lk, lij, lik, ljk, lijk
    real(dp), intent(in) :: xi, xj, xk, xij, xik, xjk, xijk
    real(dp) :: a, b, c, ak, bk, ck
    a = lij + li*lj
    b = li*xj + lj*xi + xij
    c = xi*xj
    ak = lijk + lik*lj + li*ljk
    bk = lik*xj + li*xjk + ljk*xi + lj*xik + xijk
    ck = xik*xj + xi*xjk
    value = prefactor*((lk*a + ak)*cdf_value + &
      (lk*b + a*xk + bk)*pdf_value - &
      (lk*c + b*xk + ck)*x*pdf_value + c*xk*(x*x - 1.0_dp)*pdf_value)
  end function third_derivative

  pure real(dp) function binary_greek(name, payoff, spot, rate, time, sigma, &
      dividend, d1, d2, fair_value) result(value)
    character(len=*), intent(in) :: name, payoff
    real(dp), intent(in) :: spot, rate, time, sigma, dividend
    real(dp), intent(in) :: d1, d2, fair_value
    logical :: asset, call_option
    real(dp) :: sgn, dshift, based, x, cdfx, pdfx, pref, sqrt_time, dbase
    real(dp) :: ls, lv, lt, lr, lq, lss, lsss
    real(dp) :: ds, dv, dtm, dr, dq, dss, dsv, dst, dvv, dvt, dvr
    real(dp) :: dsss, dssv, dsst, dvvv
    real(dp) :: xs, xv, xt, xr, xq, xss, xsv, xst, xvv, xvt, xvr
    real(dp) :: xsss, xssv, xsst, xvvv

    asset = payoff == 'asset_or_nothing_call' .or. payoff == 'asset_or_nothing_put'
    call_option = payoff == 'cash_or_nothing_call' .or. &
      payoff == 'digital_call' .or. payoff == 'asset_or_nothing_call'
    sgn = merge(1.0_dp, -1.0_dp, call_option)
    dshift = merge(0.5_dp, -0.5_dp, asset)
    based = merge(d1, d2, asset)
    x = sgn*based
    cdfx = normal_cdf(x)
    pdfx = normal_pdf(x)
    pref = merge(spot*exp(-dividend*time), exp(-rate*time), asset)
    sqrt_time = sqrt(time)

    ls = merge(1.0_dp/spot, 0.0_dp, asset)
    lv = 0.0_dp
    lt = merge(-dividend, -rate, asset)
    lr = merge(0.0_dp, -time, asset)
    lq = merge(-time, 0.0_dp, asset)
    lss = merge(-1.0_dp/spot**2, 0.0_dp, asset)
    lsss = merge(2.0_dp/spot**3, 0.0_dp, asset)

    dbase = based - dshift*sigma*sqrt_time
    ds = 1.0_dp/(spot*sigma*sqrt_time)
    dv = -based/sigma + 2.0_dp*dshift*sqrt_time
    dtm = (rate - dividend + dshift*sigma**2)/(sigma*sqrt_time) - based/(2.0_dp*time)
    dr = sqrt_time/sigma
    dq = -sqrt_time/sigma
    dss = -1.0_dp/(spot**2*sigma*sqrt_time)
    dsv = -1.0_dp/(spot*sigma**2*sqrt_time)
    dst = -1.0_dp/(2.0_dp*spot*sigma*time*sqrt_time)
    dvv = 2.0_dp*dbase/sigma**2
    dvt = -dtm/sigma + dshift/sqrt_time
    dvr = -sqrt_time/sigma**2
    dsss = 2.0_dp/(spot**3*sigma*sqrt_time)
    dssv = 1.0_dp/(spot**2*sigma**2*sqrt_time)
    dsst = 1.0_dp/(2.0_dp*spot**2*sigma*time*sqrt_time)
    dvvv = -6.0_dp*dbase/sigma**3

    xs = sgn*ds; xv = sgn*dv; xt = sgn*dtm; xr = sgn*dr; xq = sgn*dq
    xss = sgn*dss; xsv = sgn*dsv; xst = sgn*dst; xvv = sgn*dvv
    xvt = sgn*dvt; xvr = sgn*dvr
    xsss = sgn*dsss; xssv = sgn*dssv; xsst = sgn*dsst; xvvv = sgn*dvvv

    value = ieee_value(0.0_dp, ieee_quiet_nan)
    select case (name)
    case ('fair_value')
      value = fair_value
    case ('delta')
      value = first_derivative(pref, cdfx, pdfx, ls, xs)
    case ('vega')
      value = first_derivative(pref, cdfx, pdfx, lv, xv)
    case ('theta')
      value = -first_derivative(pref, cdfx, pdfx, lt, xt)
    case ('rho')
      value = first_derivative(pref, cdfx, pdfx, lr, xr)
    case ('epsilon')
      value = first_derivative(pref, cdfx, pdfx, lq, xq)
    case ('lambda')
      value = first_derivative(pref, cdfx, pdfx, ls, xs)*spot/fair_value
    case ('gamma')
      value = second_derivative(pref, cdfx, pdfx, x, ls, ls, lss, xs, xs, xss)
    case ('vanna')
      value = second_derivative(pref, cdfx, pdfx, x, ls, lv, 0.0_dp, xs, xv, xsv)
    case ('charm')
      value = -second_derivative(pref, cdfx, pdfx, x, ls, lt, 0.0_dp, xs, xt, xst)
    case ('vomma')
      value = second_derivative(pref, cdfx, pdfx, x, lv, lv, 0.0_dp, xv, xv, xvv)
    case ('veta')
      value = second_derivative(pref, cdfx, pdfx, x, lv, lt, 0.0_dp, xv, xt, xvt)
    case ('vera')
      value = second_derivative(pref, cdfx, pdfx, x, lr, lv, 0.0_dp, xr, xv, xvr)
    case ('speed')
      value = third_derivative(pref, cdfx, pdfx, x, ls, ls, ls, lss, lss, &
        lss, lsss, xs, xs, xs, xss, xss, xss, xsss)
    case ('zomma')
      value = third_derivative(pref, cdfx, pdfx, x, ls, ls, lv, lss, 0.0_dp, &
        0.0_dp, 0.0_dp, xs, xs, xv, xss, xsv, xsv, xssv)
    case ('color')
      value = third_derivative(pref, cdfx, pdfx, x, ls, ls, lt, lss, 0.0_dp, &
        0.0_dp, 0.0_dp, xs, xs, xt, xss, xst, xst, xsst)
    case ('ultima')
      value = third_derivative(pref, cdfx, pdfx, x, lv, lv, lv, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, xv, xv, xv, xvv, xvv, xvv, xvvv)
    end select
  end function binary_greek

  subroutine bs_geometric_asian_greeks(spot, strike, rate, time, sigma, dividend, &
      payoff, requested, result)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    character(len=*), intent(in) :: payoff
    character(len=*), intent(in) :: requested(:)
    type(greek_result), intent(out) :: result
    real(dp) :: lg, st, st3, sq3, sq3t, dg, dd, dv, dtm, dgg
    real(dp) :: shift, ad, disc, pref, dk, dvs, dv2, adv, adv2, lpv, lpv2
    integer :: i
    logical :: valid

    call initialize_result(result, requested)
    call validate_inputs(spot, strike, time, sigma, result, valid)
    if (.not. valid) return
    if (trim(payoff) /= 'call' .and. trim(payoff) /= 'put') then
      call set_error(result, greeks_unknown_payoff, &
        'geometric Asian payoff must be call or put')
      return
    end if
    lg = log(spot/strike)
    st = sqrt(time); st3 = sqrt(time/3.0_dp); sq3 = sqrt(3.0_dp); sq3t = sqrt(3.0_dp*time)
    dg = (lg + 0.5_dp*time*(rate - dividend - 0.5_dp*sigma**2))/(sigma*st3)
    dd = sq3/(spot*sigma*st)
    dv = -sq3t*(0.25_dp + lg/(sigma**2*time) + (rate - dividend)/(2.0_dp*sigma**2))
    dtm = -sq3*lg/(2.0_dp*sigma*time**1.5_dp) + &
      sq3*(rate - dividend - 0.5_dp*sigma**2)/(4.0_dp*sigma*st)
    dgg = -sq3/(spot**2*sigma*st)
    shift = sigma*st3
    ad = dg + shift
    disc = exp(-rate*time)
    pref = exp(-0.5_dp*time*(rate + dividend + sigma**2/6.0_dp))
    dk = disc*strike
    dvs = st3
    dv2 = 2.0_dp*sq3t*(lg/time + 0.5_dp*(rate - dividend))/sigma**3
    adv = dv + dvs; adv2 = dv2
    lpv = -sigma*time/6.0_dp; lpv2 = -time/6.0_dp

    do i = 1, size(requested)
      result%values(i) = geometric_value(trim(requested(i)), trim(payoff), spot, &
        strike, rate, time, sigma, dividend, dg, dd, dv, dtm, dgg, ad, disc, &
        pref, dk, adv, adv2, lpv, lpv2)
      if (.not. ieee_is_finite_local(result%values(i))) then
        call set_error(result, greeks_numerical_error, &
          'unknown geometric-Asian Greek or nonfinite result')
        return
      end if
    end do
  end subroutine bs_geometric_asian_greeks

  pure real(dp) function geometric_value(name, payoff, spot, strike, rate, time, &
      sigma, dividend, dg, dd, dv, dtm, dgg, ad, disc, pref, dk, adv, adv2, &
      lpv, lpv2) result(value)
    character(len=*), intent(in) :: name, payoff
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    real(dp), intent(in) :: dg, dd, dv, dtm, dgg, ad, disc, pref, dk
    real(dp), intent(in) :: adv, adv2, lpv, lpv2
    real(dp) :: a, b, normal_time_shift

    value = ieee_value(0.0_dp, ieee_quiet_nan)
    normal_time_shift = sigma/(2.0_dp*sqrt(3.0_dp*time))
    if (payoff == 'call') then
      select case (name)
      case ('fair_value')
        value = spot*pref*normal_cdf(ad) - dk*normal_cdf(dg)
      case ('delta')
        value = pref*(normal_cdf(ad) + spot*normal_pdf(ad)*dd) - &
          dk*normal_pdf(dg)*dd
      case ('vega')
        value = spot*pref*(lpv*normal_cdf(ad) + normal_pdf(ad)*adv) - &
          dk*normal_pdf(dg)*dv
      case ('rho')
        value = spot*pref*(-0.5_dp*time*normal_cdf(ad) + &
          normal_pdf(ad)*sqrt(3.0_dp*time)/(2.0_dp*sigma)) + &
          strike*disc*(time*normal_cdf(dg) - &
          normal_pdf(dg)*sqrt(3.0_dp*time)/(2.0_dp*sigma))
      case ('theta')
        value = -spot*pref*(-0.5_dp*(rate + dividend + sigma**2/6.0_dp) * &
          normal_cdf(ad) + normal_pdf(ad)*(dtm + normal_time_shift)) + &
          strike*disc*(-rate*normal_cdf(dg) + normal_pdf(dg)*dtm)
      case ('gamma')
        value = pref*normal_pdf(ad)*(2.0_dp*dd - spot*ad*dd**2 + spot*dgg) + &
          dk*normal_pdf(dg)*(dg*dd**2 - dgg)
      case ('vomma')
        value = spot*pref*((lpv**2 + lpv2)*normal_cdf(ad) + &
          (2.0_dp*lpv*adv - ad*adv**2 + adv2)*normal_pdf(ad)) + &
          dk*normal_pdf(dg)*(dg*dv**2 - adv2)
      end select
    else
      a = -ad
      b = -dg
      select case (name)
      case ('fair_value')
        value = -spot*pref*normal_cdf(a) + dk*normal_cdf(b)
      case ('delta')
        value = -pref*(normal_cdf(a) - spot*normal_pdf(a)*dd) - &
          dk*normal_pdf(b)*dd
      case ('vega')
        value = -spot*pref*(lpv*normal_cdf(a) - normal_pdf(a)*adv) - &
          dk*normal_pdf(b)*dv
      case ('rho')
        value = -spot*pref*(-0.5_dp*time*normal_cdf(a) - &
          normal_pdf(a)*sqrt(3.0_dp*time)/(2.0_dp*sigma)) - &
          strike*disc*(time*normal_cdf(b) + &
          normal_pdf(b)*sqrt(3.0_dp*time)/(2.0_dp*sigma))
      case ('theta')
        value = spot*pref*(-0.5_dp*(rate + dividend + sigma**2/6.0_dp) * &
          normal_cdf(a) - normal_pdf(a)*(dtm + normal_time_shift)) + &
          strike*disc*(rate*normal_cdf(b) + normal_pdf(b)*dtm)
      case ('gamma')
        value = pref*normal_pdf(a)*(2.0_dp*dd - spot*a*dd**2 + spot*dgg) - &
          dk*normal_pdf(b)*(dg*dd**2 - dgg)
      case ('vomma')
        value = -spot*pref*((lpv**2 + lpv2)*normal_cdf(a) + &
          (-2.0_dp*lpv*adv + ad*adv**2 - adv2)*normal_pdf(ad)) + &
          dk*normal_pdf(dg)*(dg*dv**2 - adv2)
      end select
    end if
  end function geometric_value

  subroutine bs_implied_volatility(option_price, spot, strike, rate, time, &
      dividend, payoff, volatility, status, iterations, start_volatility, &
      precision, max_iter)
    real(dp), intent(in) :: option_price, spot, strike, rate, time, dividend
    character(len=*), intent(in) :: payoff
    real(dp), intent(out) :: volatility
    integer, intent(out) :: status, iterations
    real(dp), intent(in), optional :: start_volatility, precision
    integer, intent(in), optional :: max_iter
    real(dp) :: tol, sigma_value, d1, d2, vega, vomma, price, error_value
    integer :: maximum, i

    status = 0; iterations = 0
    sigma_value = 0.3_dp
    if (present(start_volatility)) sigma_value = start_volatility
    tol = 1.0e-6_dp
    if (present(precision)) tol = precision
    maximum = 30
    if (present(max_iter)) maximum = max_iter
    volatility = sigma_value
    if (option_price <= bs_european_price(spot, strike, rate, time, 1.0e-12_dp, &
        dividend, payoff) .or. sigma_value <= 0.0_dp .or. tol <= 0.0_dp .or. &
        maximum <= 0) then
      status = greeks_invalid_argument
      return
    end if
    do i = 1, maximum
      d1 = (log(spot/strike) + (rate - dividend + 0.5_dp*sigma_value**2)*time) / &
        (sigma_value*sqrt(time))
      d2 = d1 - sigma_value*sqrt(time)
      price = bs_european_price(spot, strike, rate, time, sigma_value, dividend, payoff)
      error_value = price - option_price
      if (abs(error_value) < tol) then
        volatility = sigma_value
        iterations = i
        return
      end if
      vega = spot*exp(-dividend*time)*normal_pdf(d1)*sqrt(time)
      vomma = vega*d1*d2/sigma_value
      if (abs(2.0_dp*vega**2 - error_value*vomma) <= tiny(1.0_dp)) exit
      sigma_value = sigma_value - 2.0_dp*error_value*vega / &
        (2.0_dp*vega**2 - error_value*vomma)
      if (.not. ieee_is_finite_local(sigma_value) .or. sigma_value <= 0.0_dp) &
        sigma_value = max(1.0e-8_dp, 0.5_dp*volatility)
      volatility = sigma_value
    end do
    status = greeks_no_convergence
    iterations = maximum
  end subroutine bs_implied_volatility

end module greeks_black_scholes
