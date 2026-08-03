! SPDX-License-Identifier: MIT
module test_payoff_module
  use greeks, only: dp
  implicit none
contains
  pure real(dp) function custom_call(x, strike) result(value)
    real(dp), intent(in) :: x, strike
    value = max(0.0_dp, x - strike)
  end function custom_call

  pure real(dp) function custom_call_derivative(x, strike) result(value)
    real(dp), intent(in) :: x, strike
    value = merge(1.0_dp, 0.0_dp, x > strike)
  end function custom_call_derivative
end module test_payoff_module

program test_monte_carlo
  use greeks
  use test_payoff_module, only: custom_call, custom_call_derivative
  implicit none
  type(greek_result) :: result, exact, custom_result
  character(len=24) :: names(3), one(1)
  real(dp) :: exact_value

  names = [character(len=24) :: 'fair_value', 'delta', 'vega']
  call malliavin_european_greeks(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, &
    0.2_dp, 'call', names, result, paths=60000, seed=42, antithetic=.true.)
  call check(result%status == 0, 'European MC status')
  call check_close(result%values(1), 10.450583572185565_dp, 0.18_dp, &
    'European MC fair value')
  call check_close(result%values(2), 0.6368306511756191_dp, 0.035_dp, &
    'European MC delta')
  call check(result%standard_errors(1) > 0.0_dp, 'European MC standard error')

  call malliavin_european_greeks(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, &
    0.2_dp, 'custom', names, custom_result, paths=60000, seed=42, &
    antithetic=.true., payoff_fn=custom_call)
  call check_close(custom_result%values(1), result%values(1), 0.0_dp, &
    'custom European payoff callback')

  one(1) = 'fair_value'
  call bs_geometric_asian_greeks(100.0_dp, 100.0_dp, 0.03_dp, 1.0_dp, &
    0.25_dp, 0.01_dp, 'call', one, exact)
  exact_value = exact%values(1)
  call malliavin_geometric_asian_greeks(100.0_dp, 100.0_dp, 0.03_dp, &
    1.0_dp, 0.25_dp, 0.01_dp, 'call', one, result, steps=32, paths=40000, &
    seed=7, antithetic=.true.)
  call check(result%status == 0, 'geometric MC status')
  call check_close(result%values(1), exact_value, 0.25_dp, 'geometric MC value')

  names = [character(len=24) :: 'fair_value', 'delta', 'rho']
  call malliavin_asian_greeks(100.0_dp, 100.0_dp, 0.03_dp, 1.0_dp, 0.25_dp, &
    0.01_dp, 'call', names, result, steps=24, paths=12000, seed=11, &
    antithetic=.true.)
  call check(result%status == 0, 'Asian MC status')
  call check(result%values(1) > 0.0_dp, 'Asian MC positive price')
  call check(result%values(2) > 0.0_dp, 'Asian MC positive delta')

  one(1) = 'delta_d'
  call malliavin_asian_greeks(100.0_dp, 100.0_dp, 0.03_dp, 1.0_dp, 0.25_dp, &
    0.01_dp, 'custom', one, custom_result, steps=12, paths=3000, seed=13, &
    payoff_fn=custom_call, derivative_fn=custom_call_derivative)
  call check(custom_result%status == 0, 'custom Asian derivative callback')
  call check(custom_result%values(1) > 0.0_dp, 'custom Asian delta_d')

  call bs_malliavin_asian_greeks(100.0_dp, 100.0_dp, 0.03_dp, 1.0_dp, &
    0.25_dp, 0.01_dp, 'call', names, result, steps=24, paths=12000, seed=11)
  call check(result%status == 0, 'control-variate status')
  call check(result%values(1) > 0.0_dp, 'control-variate positive price')

  one(1) = 'fair_value'
  call malliavin_asian_greeks(100.0_dp, 100.0_dp, 0.03_dp, 1.0_dp, 0.25_dp, &
    0.01_dp, 'put', one, result, model='jump_diffusion', &
    jump_intensity=0.1_dp, jump_scale=0.05_dp, steps=12, paths=4000, seed=3)
  call check(result%status == 0, 'jump-diffusion status')
  call check(result%values(1) >= 0.0_dp, 'jump-diffusion price')

  print '(a)', 'test_monte_carlo: PASS'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine check

  subroutine check_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call check(abs(actual - expected) <= tolerance, label)
  end subroutine check_close
end program test_monte_carlo
