! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program test_finance
  use r4good_personal_finances
  implicit none
  real(dp), allocatable :: pv(:)
  real(dp) :: cash(10), h
  integer :: st

  call assert_close(purchasing_power(1.0_dp, 50.0_dp, 0.02_dp), 2.691588029_dp, 1.0e-8_dp, "positive purchasing power")
  call assert_close(purchasing_power(10.0_dp, 30.0_dp, -0.02_dp), 5.520708890_dp, 1.0e-8_dp, "negative purchasing power")
  call assert_close(optimal_risky_asset_allocation(0.05_dp, 0.20_dp, 0.0_dp, 2.0_dp), 0.625_dp, 1.0e-12_dp, "Merton share")
  call assert_close(optimal_risky_asset_allocation(0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp), 0.0_dp, 1.0e-12_dp, "zero Merton share")
  call assert_close(risk_adjusted_return(0.015_dp, 0.055_dp, 0.5_dp, 0.20_dp, 2.0_dp), 0.025_dp, 1.0e-12_dp, "risk adjusted return")
  call assert_close(risk_adjusted_return(0.015_dp, 0.055_dp, 0.5_dp), 0.025_dp, 1.0e-12_dp, "optimal risk adjusted shortcut")
  cash = 100000.0_dp
  call present_value_stream(cash, 0.02_dp, pv, st)
  call assert_true(st == r4gpf_success, "present value status")
  call assert_close(pv(1), 916223.671_dp, 1.0e-3_dp, "present value")
  call assert_close(utility(0.5_dp, 0.5_dp), -1.0_dp, 1.0e-12_dp, "CRRA utility")
  call assert_close(utility(1.5_dp, 1.0_dp), log(1.5_dp), 1.0e-12_dp, "log utility")
  call assert_close(inverse_utility(utility(2.0_dp, 0.2_dp), 0.2_dp), 2.0_dp, 1.0e-11_dp, "inverse utility")
  h = certainty_equivalent_return(0.03370865_dp, 0.004088745_dp, 0.25_dp)
  call assert_close(h, 0.0252899766_dp, 1.0e-8_dp, "certainty equivalent")
  call assert_close(discretionary_spending(400000.0_dp, 0.02_dp, 0.5_dp, 35.0_dp, 120.0_dp, &
    90.0_dp, 10.0_dp, discount_rate=0.03_dp), 13167.59957_dp, 1.0e-5_dp, "discretionary spending")
  print *, "test_finance: PASS"
contains
  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual - expected) > tolerance) then
      print *, "FAIL: ", trim(message), actual, expected
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      print *, "FAIL: ", trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_finance
