! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program test_mortality
  use r4good_personal_finances
  implicit none
  type(gompertz_fit) :: fit
  real(dp), allocatable :: members(:, :), joint(:)
  real(dp) :: mode
  integer :: st

  call assert_close(gompertz_survival_probability(65.0_dp, 85.0_dp, 80.0_dp, 10.0_dp), 0.2404_dp, 1.0e-3_dp, "Gompertz survival 1")
  call assert_close(gompertz_survival_probability(75.0_dp, 85.0_dp, 80.0_dp, 10.0_dp), 0.3527_dp, 1.0e-3_dp, "Gompertz survival 2")
  call assert_close(gompertz_survival_probability(100.0_dp, 100.0_dp, &
       80.0_dp, 10.0_dp, 100.0_dp), 0.0_dp, 1.0e-14_dp, &
       "truncated endpoint")
  call assert_close(life_expectancy(25.0_dp, 91.0_dp, 8.88_dp, 115.0_dp), 85.91561044_dp, 1.0e-7_dp, "life expectancy")
  call gompertz_mode_from_life_expectancy(86.0_dp, 25.0_dp, 8.88_dp, mode, st, 115.0_dp)
  call assert_true(st == r4gpf_success, "mode inversion status")
  call assert_close(mode, 91.08473352_dp, 1.0e-6_dp, "mode inversion")
  call fit_joint_gompertz([65.0_dp, 65.0_dp], [88.0_dp, 91.0_dp], [10.65_dp, 8.88_dp], 110.0_dp, fit, members, joint)
  call assert_true(fit%status == r4gpf_success, "joint fit status")
  call assert_close(fit%mode, 94.7642123_dp, 1.0e-5_dp, "joint mode")
  call assert_close(fit%dispersion, 6.1808649_dp, 1.0e-5_dp, "joint dispersion")
  call assert_close(joint(1), 1.0_dp, 1.0e-14_dp, "joint survival starts at one")
  print *, "test_mortality: PASS"
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
end program test_mortality
