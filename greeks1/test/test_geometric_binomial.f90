! SPDX-License-Identifier: MIT
program test_geometric_binomial
  use greeks
  implicit none
  type(greek_result) :: result
  character(len=24) :: names(7)
  real(dp) :: european_put

  names = [character(len=24) :: 'fair_value', 'delta', 'vega', 'rho', &
    'theta', 'gamma', 'vomma']
  call bs_geometric_asian_greeks(110.0_dp, 100.0_dp, 0.02_dp, 4.5_dp, &
    0.22_dp, 0.015_dp, 'call', names, result)
  call check_close(result%values(1), 15.01411771631934_dp, 2.0e-12_dp, &
    'geometric call value')
  call check_close(result%values(2), 0.6156218630925259_dp, 2.0e-12_dp, &
    'geometric call delta')
  call check_close(result%values(3), 32.65401327051427_dp, 2.0e-10_dp, &
    'geometric call vega')
  call check_close(result%values(4), 84.80288139196313_dp, 2.0e-10_dp, &
    'geometric call rho')
  call check_close(result%values(5), -0.6672228713032955_dp, 2.0e-12_dp, &
    'geometric call theta')
  call check_close(result%values(6), 0.010976095688866422_dp, 2.0e-13_dp, &
    'geometric call gamma')
  call check_close(result%values(7), -17.115702893382036_dp, 2.0e-10_dp, &
    'geometric call vomma')

  call bs_geometric_asian_greeks(110.0_dp, 100.0_dp, 0.02_dp, 4.5_dp, &
    0.22_dp, 0.015_dp, 'put', names, result)
  call check_close(result%values(1), 6.566091953542262_dp, 2.0e-12_dp, &
    'geometric put value')
  call check_close(result%values(2), -0.2920249031792912_dp, 2.0e-12_dp, &
    'geometric put delta')
  call check_close(result%values(3), 49.127802078347756_dp, 2.0e-10_dp, &
    'geometric put vega')
  call check_close(result%values(4), -101.8235773278148_dp, 2.0e-10_dp, &
    'geometric put rho')
  call check_close(result%values(5), -0.9892731411366834_dp, 2.0e-12_dp, &
    'geometric put theta')
  call check_close(result%values(6), 0.01097609568886642_dp, 2.0e-13_dp, &
    'geometric put gamma')
  call check_close(result%values(7), 55.04698017075032_dp, 2.0e-9_dp, &
    'geometric put vomma')

  names = [character(len=24) :: 'fair_value', 'delta', 'gamma', 'vega', &
    'theta', 'rho', 'epsilon']
  call binomial_american_greeks(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, &
    0.2_dp, 0.0_dp, 'call', names, result, steps=300)
  call check_close(result%values(1), 10.450583572185565_dp, 2.0e-6_dp, &
    'American nondividend call equals European')

  european_put = bs_european_price(100.0_dp, 105.0_dp, 0.05_dp, 1.0_dp, &
    0.2_dp, 0.0_dp, 'put')
  call binomial_american_greeks(100.0_dp, 105.0_dp, 0.05_dp, 1.0_dp, &
    0.2_dp, 0.0_dp, 'put', names, result, steps=500)
  call check(result%values(1) >= european_put, 'American put premium')
  call check(result%values(2) < 0.0_dp, 'American put delta')

  print '(a)', 'test_geometric_binomial: PASS'
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
end program test_geometric_binomial
