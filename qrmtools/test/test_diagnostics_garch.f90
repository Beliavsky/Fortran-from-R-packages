! SPDX-License-Identifier: GPL-3.0-or-later
program test_diagnostics_garch
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use qrmtools, only : dp, mahalanobis_squared, maha2_test, mardia_test, &
    test_result, loglik_garch_11, fit_garch_11, tail_index_garch_11, &
    garch_result
  implicit none

  real(dp) :: x(8,2)
  real(dp), allocatable :: distances(:)
  type(test_result) :: test
  type(garch_result) :: fit
  real(dp) :: series(5)
  real(dp) :: variance0
  real(dp) :: tail

  x = reshape([-1.0_dp,-0.7_dp,-0.2_dp,0.0_dp,0.3_dp,0.7_dp,1.0_dp,1.3_dp, &
                0.2_dp,-0.4_dp,0.5_dp,-0.1_dp,0.8_dp,-0.6_dp,0.4_dp,1.1_dp],[8,2])
  distances = mahalanobis_squared(x)
  call assert_true(all(distances>=0.0_dp))
  call assert_close(sum(distances),14.0_dp,3.0e-11_dp)

  test = maha2_test(x)
  call assert_true(test%ok)
  call assert_true(test%p_value>=0.0_dp .and. test%p_value<=1.0_dp)
  test = mardia_test(x)
  call assert_true(test%ok)
  test = mardia_test(x,.true.)
  call assert_true(test%ok)

  series = [0.1_dp,-0.2_dp,0.05_dp,0.12_dp,-0.08_dp]
  variance0 = sum(series*series)/real(size(series),dp)
  call assert_close(loglik_garch_11([1.0_dp,1.0_dp],series,variance0,1.0_dp), &
    3.277789384185003_dp,3.0e-11_dp)

  fit = fit_garch_11(series,initial=[1.0_dp,1.0_dp],max_iterations=1000)
  call assert_true(fit%ok)
  call assert_true(all(fit%sigma>0.0_dp))
  call assert_true(fit%coefficients(1)>0.0_dp)
  call assert_true(fit%coefficients(2)>=0.0_dp)
  call assert_true(fit%coefficients(3)>=0.0_dp)
  call assert_true(sum(fit%coefficients(2:3))<1.0_dp)
  call assert_true(ieee_is_finite(fit%log_likelihood))

  tail = tail_index_garch_11([-1.5_dp,-1.0_dp,-0.5_dp,0.0_dp,0.5_dp, &
    1.0_dp,1.5_dp],0.08_dp,0.85_dp)
  call assert_true(ieee_is_finite(tail))
  call assert_true(tail>0.0_dp)

  print '(a)', 'test_diagnostics_garch: PASS'

contains

  subroutine assert_close(actual,expected,tolerance)
    real(dp), intent(in) :: actual,expected,tolerance
    if(abs(actual-expected)>tolerance*max(1.0_dp,abs(expected))) then
      print *, 'mismatch:',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if(.not.condition) error stop 1
  end subroutine assert_true

end program test_diagnostics_garch
