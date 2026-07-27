! SPDX-License-Identifier: MIT
program test_anchored
  use epo, only : dp, epo_from_covariance, epo_result
  implicit none

  real(dp) :: anchor(3), covariance(3,3), signal(3), simple_weights(3)
  type(epo_result) :: fit

  covariance = reshape([ &
     0.040_dp,  0.006_dp, -0.004_dp, &
     0.006_dp,  0.090_dp,  0.012_dp, &
    -0.004_dp,  0.012_dp,  0.025_dp  &
  ], [3,3])
  signal = [0.08_dp, 0.11_dp, 0.05_dp]
  anchor = [0.4_dp, 0.3_dp, 0.3_dp]

  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'anchored', &
    0.35_dp, anchor=anchor, normalize=.false., endogenous=.true.)
  call assert_true(fit%ok)
  call assert_close(fit%gamma, 0.2251169670523295_dp, 3.0e-14_dp)
  call assert_vector(fit%weights, [ &
    0.4339951239928575_dp, &
    0.2329608076862335_dp, &
    0.3701037780651807_dp  &
  ], 3.0e-14_dp)

  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'anchored', &
    0.35_dp, anchor=anchor, normalize=.true., endogenous=.true.)
  call assert_vector(fit%weights, [ &
    0.4184861487868199_dp, &
    0.2246358676335804_dp, &
    0.3568779835795998_dp  &
  ], 4.0e-14_dp)

  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'anchored', &
    0.35_dp, anchor=anchor, normalize=.true., endogenous=.false.)
  call assert_vector(fit%weights, [ &
    0.4172644225850857_dp, &
    0.2178420961482155_dp, &
    0.3648934812666987_dp  &
  ], 4.0e-14_dp)

  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'simple', &
    0.0_dp, normalize=.true.)
  simple_weights = fit%weights
  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'anchored', &
    0.0_dp, anchor=anchor, normalize=.true., endogenous=.true.)
  call assert_vector(fit%weights, simple_weights, 2.0e-14_dp)

  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'anchored', &
    1.0_dp, anchor=anchor, normalize=.true., endogenous=.true.)
  call assert_vector(fit%weights, anchor, 2.0e-14_dp)

  print '(a)', 'test_anchored: PASS'

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

end program test_anchored
