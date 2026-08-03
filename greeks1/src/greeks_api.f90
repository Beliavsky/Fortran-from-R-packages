! SPDX-License-Identifier: MIT
module greeks_api
  use greeks_kinds, only: dp
  use greeks_types, only: greek_result, initialize_result, set_error
  use greeks_types, only: greeks_invalid_argument, greeks_no_convergence
  use greeks_black_scholes, only: bs_european_greeks, bs_geometric_asian_greeks
  use greeks_black_scholes, only: bs_implied_volatility
  use greeks_binomial, only: binomial_american_greeks
  use greeks_monte_carlo, only: malliavin_asian_greeks
  use greeks_monte_carlo, only: bs_malliavin_asian_greeks
  implicit none
  private
  public :: option_greeks, implied_volatility

contains

  subroutine option_greeks(spot, strike, rate, time, sigma, dividend, model, &
      option_type, payoff, requested, result, steps, paths, seed, antithetic)
    real(dp), intent(in) :: spot, strike, rate, time, sigma, dividend
    character(len=*), intent(in) :: model, option_type, payoff
    character(len=*), intent(in) :: requested(:)
    type(greek_result), intent(out) :: result
    integer, intent(in), optional :: steps, paths, seed
    logical, intent(in), optional :: antithetic
    character(len=:), allocatable :: type_key, model_key
    integer :: nsteps, npaths, seed_value
    logical :: anti, basic_asian
    integer :: i

    type_key = lower_ascii(trim(option_type))
    model_key = lower_ascii(trim(model))
    if (type_key == 'american') then
      nsteps = 1000
    else
      nsteps = max(1, nint(time*252.0_dp))
    end if
    if (present(steps)) nsteps = steps
    npaths = 10000; if (present(paths)) npaths = paths
    seed_value = 1; if (present(seed)) seed_value = seed
    anti = .true.; if (present(antithetic)) anti = antithetic

    if ((type_key == 'european' .or. type_key == 'digital') .and. &
        model_key == 'black_scholes') then
      call bs_european_greeks(spot, strike, rate, time, sigma, dividend, payoff, &
        requested, result)
    else if (type_key == 'american' .and. model_key == 'black_scholes') then
      call binomial_american_greeks(spot, strike, rate, time, sigma, dividend, &
        payoff, requested, result, nsteps)
    else if (type_key == 'geometric asian' .and. model_key == 'black_scholes') then
      call bs_geometric_asian_greeks(spot, strike, rate, time, sigma, dividend, &
        payoff, requested, result)
    else if (type_key == 'asian' .and. model_key == 'black_scholes') then
      basic_asian = .true.
      do i = 1, size(requested)
        if (trim(requested(i)) /= 'fair_value' .and. trim(requested(i)) /= 'delta' .and. &
            trim(requested(i)) /= 'rho' .and. trim(requested(i)) /= 'vega') then
          basic_asian = .false.
        end if
      end do
      if (basic_asian) then
        call bs_malliavin_asian_greeks(spot, strike, rate, time, sigma, dividend, &
          payoff, requested, result, nsteps, npaths, seed_value)
      else
        call malliavin_asian_greeks(spot, strike, rate, time, sigma, dividend, &
          payoff, requested, result, model_key, steps=nsteps, paths=npaths, &
          seed=seed_value, antithetic=anti)
      end if
    else if (type_key == 'asian') then
      call malliavin_asian_greeks(spot, strike, rate, time, sigma, dividend, &
        payoff, requested, result, model_key, steps=nsteps, paths=npaths, &
        seed=seed_value, antithetic=anti)
    else
      call initialize_result(result, requested)
      call set_error(result, greeks_invalid_argument, &
        'unknown model/option-type combination')
    end if
  end subroutine option_greeks

  subroutine implied_volatility(option_price, spot, strike, rate, time, dividend, &
      model, option_type, payoff, volatility, status, iterations, &
      start_volatility, precision, max_iter, steps, paths, seed)
    real(dp), intent(in) :: option_price, spot, strike, rate, time, dividend
    character(len=*), intent(in) :: model, option_type, payoff
    real(dp), intent(out) :: volatility
    integer, intent(out) :: status, iterations
    real(dp), intent(in), optional :: start_volatility, precision
    integer, intent(in), optional :: max_iter, steps, paths, seed
    real(dp) :: sigma_value, tol, price_value, vega_value
    integer :: maximum, nsteps, npaths, seed_value, i
    type(greek_result), allocatable :: result
    character(len=24) :: requested(2)

    if (lower_ascii(trim(option_type)) == 'european' .and. &
        (trim(payoff) == 'call' .or. trim(payoff) == 'put')) then
      call bs_implied_volatility(option_price, spot, strike, rate, time, dividend, &
        payoff, volatility, status, iterations, start_volatility, precision, max_iter)
      return
    end if

    sigma_value = 0.3_dp; if (present(start_volatility)) sigma_value = start_volatility
    tol = 1.0e-6_dp; if (present(precision)) tol = precision
    maximum = 30; if (present(max_iter)) maximum = max_iter
    nsteps = max(1, nint(time*252.0_dp)); if (present(steps)) nsteps = steps
    npaths = 10000; if (present(paths)) npaths = paths
    seed_value = 1; if (present(seed)) seed_value = seed
    requested = [character(len=24) :: 'fair_value', 'vega']
    allocate(result)
    status = 0; iterations = 0; volatility = sigma_value
    if (option_price < 0.0_dp .or. sigma_value <= 0.0_dp .or. tol <= 0.0_dp .or. &
        maximum <= 0) then
      status = greeks_invalid_argument
      return
    end if

    do i = 1, maximum
      call option_greeks(spot, strike, rate, time, sigma_value, dividend, model, &
        option_type, payoff, requested, result, nsteps, npaths, seed_value, .true.)
      if (result%status /= 0) then
        status = result%status
        iterations = i
        return
      end if
      price_value = result%values(1)
      vega_value = result%values(2)
      if (abs(price_value - option_price) < tol) then
        volatility = sigma_value
        iterations = i
        return
      end if
      if (abs(vega_value) <= sqrt(tiny(1.0_dp))) then
        status = greeks_no_convergence
        iterations = i
        return
      end if
      sigma_value = max(1.0e-6_dp, sigma_value - &
        (price_value - option_price)/vega_value)
      volatility = sigma_value
    end do
    status = greeks_no_convergence
    iterations = maximum
  end subroutine implied_volatility

  pure function lower_ascii(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower(i:i) = achar(code + iachar('a') - iachar('A'))
      else
        lower(i:i) = text(i:i)
      end if
    end do
  end function lower_ascii

end module greeks_api
