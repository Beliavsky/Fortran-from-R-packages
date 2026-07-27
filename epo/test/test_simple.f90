! SPDX-License-Identifier: MIT
program test_simple
  use epo, only : dp, epo_from_covariance, epo_result
  implicit none

  real(dp) :: covariance(3,3), signal(3)
  type(epo_result) :: fit

  covariance = reshape([ &
     0.040_dp,  0.006_dp, -0.004_dp, &
     0.006_dp,  0.090_dp,  0.012_dp, &
    -0.004_dp,  0.012_dp,  0.025_dp  &
  ], [3,3])
  signal = [0.08_dp, 0.11_dp, 0.05_dp]

  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'simple', &
    0.35_dp, normalize=.false.)
  call assert_true(fit%ok)
  call assert_vector(fit%weights, [ &
    0.6765366074622499_dp, &
    0.3229476991356333_dp, &
    0.6362667917124231_dp  &
  ], 2.0e-14_dp)

  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'SIMPLE', &
    0.35_dp, normalize=.true.)
  call assert_true(fit%ok)
  call assert_vector(fit%weights, [ &
    0.4135938579904344_dp, &
    0.1974308312977634_dp, &
    0.3889753107118023_dp  &
  ], 2.0e-14_dp)
  call assert_close(sum(fit%weights), 1.0_dp, 2.0e-15_dp)

  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'simple', &
    0.0_dp, normalize=.true.)
  call assert_vector(fit%weights, [ &
    0.4284405071024897_dp, &
    0.1710707194134719_dp, &
    0.4004887734840385_dp  &
  ], 2.0e-14_dp)

  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'simple', &
    1.0_dp, normalize=.true.)
  call assert_vector(fit%weights, [ &
    0.3829787234042554_dp, &
    0.2340425531914894_dp, &
    0.3829787234042554_dp  &
  ], 2.0e-14_dp)

  print '(a)', 'test_simple: PASS'

contains

  subroutine assert_vector(actual, expected, tolerance)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    integer :: i

    if (size(actual) /= size(expected)) error stop 1
    do i = 1, size(actual)
      call assert_close(actual(i), expected(i), tolerance)
    end do
  end subroutine assert_vector

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

end program test_simple
