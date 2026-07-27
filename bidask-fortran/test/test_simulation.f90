! SPDX-License-Identifier: MIT
program test_simulation
  use iso_fortran_env, only: int64
  use bidask, only: dp, ohlc_data, simulate_ohlc, edge_estimate
  implicit none
  type(ohlc_data) :: a, b, empty, large
  real(dp) :: estimate
  integer :: i

  a = simulate_ohlc(25, 20, probability=0.8_dp, spread=0.01_dp, &
    volatility=0.03_dp, overnight=0.005_dp, drift=0.001_dp, seed=12345_int64)
  b = simulate_ohlc(25, 20, probability=0.8_dp, spread=0.01_dp, &
    volatility=0.03_dp, overnight=0.005_dp, drift=0.001_dp, seed=12345_int64)
  call assert_true(maxval(abs(a%open - b%open)) <= 0.0_dp)
  call assert_true(maxval(abs(a%high - b%high)) <= 0.0_dp)
  do i = 1, a%size()
    call assert_true(a%high(i) >= max(a%open(i), a%close(i)))
    call assert_true(a%low(i) <= min(a%open(i), a%close(i)))
  end do

  empty = simulate_ohlc(10, 5, probability=0.0_dp, spread=0.01_dp, seed=77_int64)
  call assert_true(maxval(abs(empty%open - empty%close)) <= 0.0_dp)
  call assert_true(maxval(abs(empty%high - empty%low)) <= 0.0_dp)
  call assert_true(maxval(abs(empty%close - empty%close(1))) <= 0.0_dp)

  large = simulate_ohlc(500, 50, probability=1.0_dp, spread=0.01_dp, &
    volatility=0.03_dp, seed=998877_int64)
  estimate = edge_estimate(large%open, large%high, large%low, large%close)
  call assert_true(estimate > 0.0_dp .and. estimate < 0.03_dp)
  print '(a)', 'test_simulation: PASS'

contains

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_simulation
