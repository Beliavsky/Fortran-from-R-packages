! SPDX-License-Identifier: MIT
program test_implied_helpers
  use greeks
  implicit none
  real(dp) :: price, volatility
  integer :: status, iterations
  real(dp), allocatable :: w(:, :), values(:), matrix(:, :)
  real(dp) :: increments(6)
  type(greek_result) :: result
  character(len=24) :: names(2)

  price = bs_european_price(120.0_dp, 100.0_dp, 0.02_dp, 2.0_dp, 0.35_dp, &
    0.01_dp, 'call')
  call bs_implied_volatility(price, 120.0_dp, 100.0_dp, 0.02_dp, 2.0_dp, &
    0.01_dp, 'call', volatility, status, iterations, 0.2_dp, 1.0e-12_dp, 40)
  call check(status == 0, 'BS implied status')
  call check_close(volatility, 0.35_dp, 2.0e-12_dp, 'BS implied recovery')

  increments = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
  call make_bm(increments, 2, 3, w)
  call check_close(w(1, 4), 9.0_dp, 0.0_dp, 'Brownian layout row 1')
  call check_close(w(2, 4), 12.0_dp, 0.0_dp, 'Brownian layout row 2')
  allocate(matrix(1, 3))
  matrix(1, :) = [1.0_dp, 2.0_dp, 3.0_dp]
  values = calc_i(matrix, 0.5_dp)
  call check_close(values(1), 2.0_dp, 0.0_dp, 'trapezoid integral')
  values = calc_i_1(matrix, 0.5_dp)
  call check_close(values(1), 1.25_dp, 0.0_dp, 'first time moment')

  names = [character(len=24) :: 'fair_value', 'delta']
  call option_greeks(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.2_dp, 0.0_dp, &
    'Black_Scholes', 'European', 'call', names, result)
  call check(result%status == 0, 'dispatcher status')
  call check_close(result%values(1), 10.450583572185565_dp, 2.0e-12_dp, &
    'dispatcher value')

  print '(a)', 'test_implied_helpers: PASS'
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
end program test_implied_helpers
