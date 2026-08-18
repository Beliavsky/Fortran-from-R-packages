! SPDX-License-Identifier: MIT
program test_estimation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
  use ewens, only : dp, ewens_mle, ewens_mle_nk, ewens_score
  implicit none
  integer :: labels(10), one(7), unique_labels(6), singleton(1)
  real(dp) :: theta

  labels = [1, 1, 1, 1, 2, 2, 3, 3, 4, 4]
  theta = ewens_mle(labels)
  call assert_close(theta, 1.95635790111585_dp, 2.0e-12_dp)
  call assert_close(ewens_score(theta, 10, 4), 0.0_dp, 2.0e-12_dp)

  one = 9
  call assert_close(ewens_mle(one), 0.0_dp, 0.0_dp)

  unique_labels = [1, 2, 3, 4, 5, 6]
  theta = ewens_mle(unique_labels)
  if (ieee_is_finite(theta) .or. theta <= 0.0_dp) error stop 2

  singleton = 1
  if (.not. ieee_is_nan(ewens_mle(singleton))) error stop 3
  if (.not. ieee_is_nan(ewens_mle_nk(5, 0))) error stop 4

  print '(a)', 'test_estimation: PASS'
contains
  subroutine assert_close(x, y, tol)
    real(dp), intent(in) :: x, y, tol
    if (abs(x - y) > tol) error stop 1
  end subroutine assert_close
end program test_estimation
