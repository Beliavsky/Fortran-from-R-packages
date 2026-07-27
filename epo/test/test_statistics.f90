! SPDX-License-Identifier: MIT
program test_statistics
  use epo, only : covariance_to_correlation, dp, sample_covariance
  implicit none

  real(dp) :: covariance(2,2), correlation(2,2), means(2)
  real(dp) :: x(4,2)
  logical :: ok

  x = reshape([ &
    1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
    2.0_dp, 4.0_dp, 4.0_dp, 8.0_dp  &
  ], [4,2])

  call sample_covariance(x, means, covariance, ok)
  call assert_true(ok)
  call assert_close(means(1), 2.5_dp, 1.0e-14_dp)
  call assert_close(means(2), 4.5_dp, 1.0e-14_dp)
  call assert_close(covariance(1,1), 5.0_dp / 3.0_dp, 1.0e-14_dp)
  call assert_close(covariance(1,2), 3.0_dp, 1.0e-14_dp)
  call assert_close(covariance(2,1), 3.0_dp, 1.0e-14_dp)
  call assert_close(covariance(2,2), 19.0_dp / 3.0_dp, 1.0e-14_dp)

  call covariance_to_correlation(covariance, correlation, ok)
  call assert_true(ok)
  call assert_close(correlation(1,1), 1.0_dp, 0.0_dp)
  call assert_close(correlation(2,2), 1.0_dp, 0.0_dp)
  call assert_close(correlation(1,2), 0.9233805168766387_dp, 1.0e-14_dp)
  call assert_close(correlation(2,1), correlation(1,2), 0.0_dp)

  print '(a)', 'test_statistics: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance

    if (abs(actual - expected) > tolerance) then
      print '(a,3es24.16)', 'mismatch: ', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition

    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_statistics
