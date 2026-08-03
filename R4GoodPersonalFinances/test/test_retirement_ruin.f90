! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program test_retirement_ruin
  use r4good_personal_finances
  implicit none
  call assert_close(upper_incomplete_gamma(2.0_dp, 3.0_dp), 0.1991482735_dp, 1.0e-9_dp, "upper gamma positive")
  call assert_close(upper_incomplete_gamma(-0.5_dp, 1.0_dp), 0.1781477118_dp, 1.0e-9_dp, "upper gamma negative")
  call assert_close(gompertz_annuity_factor(0.01_dp, 65.0_dp, 88.0_dp, 10.0_dp), 17.88846_dp, 1.0e-4_dp, "annuity factor")
  call assert_close(retirement_ruin_probability(0.02_dp, 0.20_dp, 65.0_dp, 88.0_dp, 10.0_dp, 0.05_dp), &
    0.3362243_dp, 1.0e-6_dp, "ruin 5 percent")
  call assert_close(retirement_ruin_probability(0.02_dp, 0.20_dp, 65.0_dp, 88.0_dp, 10.0_dp, 0.04_dp), &
    0.2216641_dp, 1.0e-6_dp, "ruin 4 percent")
  call assert_close(retirement_ruin_probability(0.02_dp, 0.15_dp, 65.0_dp, 88.0_dp, 10.0_dp, 0.04_dp), &
    0.1744054_dp, 1.0e-6_dp, "ruin lower volatility")
  print *, "test_retirement_ruin: PASS"
contains
  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual - expected) > tolerance) then
      print *, "FAIL: ", trim(message), actual, expected
      error stop 1
    end if
  end subroutine assert_close
end program test_retirement_ruin
