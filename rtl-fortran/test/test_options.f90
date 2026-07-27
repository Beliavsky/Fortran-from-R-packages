! SPDX-License-Identifier: MIT
program test_options
  use rtl, only: dp, option_result, spread_option_result, crr_tree_result
  use rtl, only: gbs_option, crr_option, crr_euro_tree, spread_option, barrier_spread_option
  implicit none

  type(option_result) :: call_result, put_result, tree_price
  type(spread_option_result) :: spread_call, spread_put, barrier_result
  type(crr_tree_result) :: tree

  call_result = gbs_option(100.0_dp, 100.0_dp, 1.0_dp, 0.05_dp, 0.02_dp, 0.2_dp, "call")
  put_result = gbs_option(100.0_dp, 100.0_dp, 1.0_dp, 0.05_dp, 0.02_dp, 0.2_dp, "put")
  call assert_close(call_result%price, 8.652528553942709_dp, 1.0e-12_dp)
  call assert_close(call_result%delta, 0.562139997789784_dp, 1.0e-12_dp)
  call assert_close(call_result%gamma, 0.018974281789763_dp, 1.0e-13_dp)
  call assert_close(put_result%price, 6.730917649163303_dp, 1.0e-12_dp)
  call assert_close(call_result%price - put_result%price, &
    100.0_dp * exp(-0.03_dp) - 100.0_dp * exp(-0.05_dp), 1.0e-12_dp)

  tree_price = crr_option(100.0_dp, 100.0_dp, 0.2_dp, 0.05_dp, 0.05_dp, &
    1.0_dp, 1000, "call", "european")
  call assert_close(tree_price%price, 10.4506_dp, 3.0e-3_dp)
  tree = crr_euro_tree(100.0_dp, 100.0_dp, 0.2_dp, 0.05_dp, 1.0_dp, 5, "call")
  call assert_true(tree%status%ok)
  call assert_close(tree%price, tree%option(0, 0), 0.0_dp)
  call assert_close(tree%asset(0, 0), 100.0_dp, 0.0_dp)

  spread_call = spread_option(100.0_dp, 110.0_dp, 5.0_dp, 0.2_dp, 0.25_dp, &
    0.5_dp, 1.0_dp, 0.05_dp, "call")
  spread_put = spread_option(100.0_dp, 110.0_dp, 5.0_dp, 0.2_dp, 0.25_dp, &
    0.5_dp, 1.0_dp, 0.05_dp, "put")
  call assert_close(spread_call%price, 7.529346008053384_dp, 2.0e-12_dp)
  call assert_close(spread_put%price, 2.773198885549821_dp, 2.0e-12_dp)
  call assert_close(spread_call%price - spread_put%price, exp(-0.05_dp) * 5.0_dp, 2.0e-10_dp)

  barrier_result = barrier_spread_option(-12.0_dp, -5.0_dp, 5.0_dp, 9.0_dp, &
    0.4_dp, 0.4_dp, 0.8_dp, 1.0e-12_dp, 0.045_dp, "call", "uo", "terminal")
  call assert_close(barrier_result%price, 2.0_dp, 1.0e-12_dp)
  call assert_true(.not. barrier_result%monitoring_used)
  barrier_result = barrier_spread_option(-15.0_dp, -5.0_dp, 5.0_dp, 9.0_dp, &
    0.4_dp, 0.4_dp, 0.8_dp, 0.25_dp, 0.045_dp, "call", "uo")
  call assert_close(barrier_result%price, 0.0_dp, 0.0_dp)

  print '(a)', 'test_options: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print '(a,3es24.15)', 'mismatch: ', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_options
