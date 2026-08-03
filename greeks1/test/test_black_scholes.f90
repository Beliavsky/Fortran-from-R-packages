! SPDX-License-Identifier: MIT
program test_black_scholes
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use greeks
  implicit none
  type(greek_result) :: result
  character(len=24) :: names(6)
  real(dp) :: call_value, put_value, parity

  names = [character(len=24) :: 'fair_value', 'delta', 'gamma', 'vega', &
    'rho', 'theta']
  call bs_european_greeks(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.2_dp, &
    0.0_dp, 'call', names, result)
  call check(result%status == 0, 'call status')
  call check_close(result%values(1), 10.450583572185565_dp, 2.0e-12_dp, &
    'call fair value')
  call check_close(result%values(2), 0.6368306511756191_dp, 2.0e-12_dp, &
    'call delta')
  call check_close(result%values(3), 0.018762017345846895_dp, 2.0e-12_dp, &
    'call gamma')
  call check_close(result%values(4), 37.52403469169379_dp, 2.0e-11_dp, &
    'call vega')
  call_value = result%values(1)

  call bs_european_greeks(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.2_dp, &
    0.0_dp, 'put', names, result)
  call check_close(result%values(1), 5.573526022256971_dp, 2.0e-12_dp, &
    'put fair value')
  put_value = result%values(1)
  parity = call_value - put_value - (100.0_dp - 100.0_dp*exp(-0.05_dp))
  call check_close(parity, 0.0_dp, 2.0e-12_dp, 'put-call parity')

  names = [character(len=24) :: 'fair_value', 'delta', 'gamma', 'vega', &
    'speed', 'ultima']
  call bs_european_greeks(100.0_dp, 105.0_dp, 0.03_dp, 1.4_dp, 0.27_dp, &
    0.01_dp, 'cash_or_nothing_call', names, result)
  call check(result%status == 0, 'binary status')
  call check(all(ieee_is_finite(result%values)), 'binary finite values')
  call check(result%values(1) > 0.0_dp .and. result%values(1) < 1.0_dp, &
    'binary price range')

  print '(a)', 'test_black_scholes: PASS'
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
end program test_black_scholes
