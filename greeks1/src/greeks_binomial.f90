! SPDX-License-Identifier: MIT
module greeks_binomial
  use greeks_kinds, only: dp
  use greeks_types, only: greek_result, initialize_result, set_error
  use greeks_types, only: greeks_invalid_argument, greeks_unknown_payoff
  use greeks_black_scholes, only: bs_european_price
  implicit none
  private
  public :: binomial_american_greeks, binomial_values

contains

  subroutine binomial_values(spot, strike, rate, time, sigma, dividend, payoff, &
      steps, american_value, european_value, status)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    character(len=*), intent(in) :: payoff
    integer, intent(in) :: steps
    real(dp), intent(out) :: american_value, european_value
    integer, intent(out) :: status
    real(dp), allocatable :: price(:), av(:), ev(:)
    real(dp) :: dt, up, down, p, q, discount_terminal
    integer :: i, j, index_value

    status = 0
    american_value = 0.0_dp
    european_value = 0.0_dp
    if (spot <= 0.0_dp .or. strike <= 0.0_dp .or. time <= 0.0_dp .or. &
        sigma <= 0.0_dp .or. steps < 1) then
      status = greeks_invalid_argument
      return
    end if
    if (trim(payoff) /= 'call' .and. trim(payoff) /= 'put') then
      status = greeks_unknown_payoff
      return
    end if

    dt = time/real(steps, dp)
    up = exp(sigma*sqrt(dt))
    down = 1.0_dp/up
    p = (exp((rate - dividend)*dt) - down)/(up - down)
    q = 1.0_dp - p
    if (p < 0.0_dp .or. p > 1.0_dp) then
      status = greeks_invalid_argument
      return
    end if
    discount_terminal = exp(-rate*time)
    allocate(price(2*steps + 1), av(steps + 1), ev(steps + 1))
    do j = 0, 2*steps
      price(j + 1) = spot*exp(-sigma*real(steps - j, dp)*sqrt(dt))
    end do
    do i = 0, steps
      index_value = 2*steps - 2*i + 1
      ev(i + 1) = discount_terminal*payoff_value_local(price(index_value), strike, payoff)
    end do
    av = ev
    do j = steps - 1, 0, -1
      do i = 0, j
        ev(i + 1) = p*ev(i + 1) + q*ev(i + 2)
        av(i + 1) = p*av(i + 1) + q*av(i + 2)
        index_value = 2*steps - 2*i + j - steps + 1
        av(i + 1) = max(exp(-rate*real(j, dp)*dt)* &
          payoff_value_local(price(index_value), strike, payoff), av(i + 1))
      end do
    end do
    american_value = av(1)
    european_value = ev(1)
  end subroutine binomial_values

  pure real(dp) function payoff_value_local(x, strike, payoff) result(value)
    real(dp), intent(in) :: x, strike
    character(len=*), intent(in) :: payoff
    if (trim(payoff) == 'call') then
      value = max(0.0_dp, x - strike)
    else
      value = max(0.0_dp, strike - x)
    end if
  end function payoff_value_local

  real(dp) function corrected_price(spot, strike, rate, time, sigma, dividend, &
      payoff, steps, status) result(value)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    character(len=*), intent(in) :: payoff
    integer, intent(in) :: steps
    integer, intent(out) :: status
    real(dp) :: american_value, european_tree, european_exact

    call binomial_values(spot, strike, rate, time, sigma, dividend, payoff, &
      steps, american_value, european_tree, status)
    if (status /= 0) then
      value = 0.0_dp
      return
    end if
    european_exact = bs_european_price(spot, strike, rate, time, sigma, &
      dividend, payoff)
    value = american_value + european_exact - european_tree
  end function corrected_price

  subroutine binomial_american_greeks(spot, strike, rate, time, sigma, dividend, &
      payoff, requested, result, steps, eps)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    character(len=*), intent(in) :: payoff
    character(len=*), intent(in) :: requested(:)
    type(greek_result), intent(out) :: result
    integer, intent(in), optional :: steps
    real(dp), intent(in), optional :: eps
    integer :: nsteps, i, status
    real(dp) :: h, hg, v0, vu, vd

    call initialize_result(result, requested)
    nsteps = 1000
    if (present(steps)) nsteps = steps
    h = 1.0e-5_dp
    if (present(eps)) h = eps
    if (h <= 0.0_dp) then
      call set_error(result, greeks_invalid_argument, 'eps must be positive')
      return
    end if
    v0 = corrected_price(spot, strike, rate, time, sigma, dividend, payoff, &
      nsteps, status)
    if (status /= 0) then
      call set_error(result, status, 'invalid American-option inputs')
      return
    end if

    do i = 1, size(requested)
      select case (trim(requested(i)))
      case ('fair_value')
        result%values(i) = v0
      case ('delta')
        vu = corrected_price(spot + h, strike, rate, time, sigma, dividend, &
          payoff, nsteps, status)
        vd = corrected_price(spot - h, strike, rate, time, sigma, dividend, &
          payoff, nsteps, status)
        result%values(i) = (vu - vd)/(2.0_dp*h)
      case ('gamma')
        hg = spot/50.0_dp
        vu = corrected_price(spot + hg, strike, rate, time, sigma, dividend, &
          payoff, nsteps, status)
        vd = corrected_price(spot - hg, strike, rate, time, sigma, dividend, &
          payoff, nsteps, status)
        result%values(i) = (vu - 2.0_dp*v0 + vd)/hg**2
      case ('vega')
        vu = corrected_price(spot, strike, rate, time, sigma + h, dividend, &
          payoff, nsteps, status)
        vd = corrected_price(spot, strike, rate, time, sigma - h, dividend, &
          payoff, nsteps, status)
        result%values(i) = (vu - vd)/(2.0_dp*h)
      case ('theta')
        vu = corrected_price(spot, strike, rate, time + h, sigma, dividend, &
          payoff, nsteps, status)
        vd = corrected_price(spot, strike, rate, time - h, sigma, dividend, &
          payoff, nsteps, status)
        result%values(i) = -(vu - vd)/(2.0_dp*h)
      case ('rho')
        vu = corrected_price(spot, strike, rate + h, time, sigma, dividend, &
          payoff, nsteps, status)
        vd = corrected_price(spot, strike, rate - h, time, sigma, dividend, &
          payoff, nsteps, status)
        result%values(i) = (vu - vd)/(2.0_dp*h)
      case ('epsilon')
        vu = corrected_price(spot, strike, rate, time, sigma, dividend + h, &
          payoff, nsteps, status)
        vd = corrected_price(spot, strike, rate, time, sigma, dividend - h, &
          payoff, nsteps, status)
        result%values(i) = (vu - vd)/(2.0_dp*h)
      case default
        call set_error(result, greeks_invalid_argument, 'unsupported American Greek')
        return
      end select
    end do
  end subroutine binomial_american_greeks

end module greeks_binomial
