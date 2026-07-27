! SPDX-License-Identifier: MIT
program test_validation
  use epo, only : dp, epo_from_covariance, epo_invalid_input, &
    epo_normalization_failure, epo_result, epo_singular_matrix
  implicit none

  real(dp) :: anchor(2), covariance(2,2), signal(2), singular(2,2)
  type(epo_result) :: fit

  covariance = reshape([0.04_dp, 0.01_dp, 0.01_dp, 0.09_dp], [2,2])
  signal = [0.08_dp, 0.05_dp]
  anchor = [0.5_dp, 0.5_dp]

  fit = epo_from_covariance(covariance, signal, 2.0_dp, 'unknown', 0.5_dp)
  call assert_status(fit, epo_invalid_input)

  fit = epo_from_covariance(covariance, signal, 2.0_dp, 'simple', -0.1_dp)
  call assert_status(fit, epo_invalid_input)

  fit = epo_from_covariance(covariance, signal, 0.0_dp, 'simple', 0.5_dp)
  call assert_status(fit, epo_invalid_input)

  fit = epo_from_covariance(covariance, signal, 2.0_dp, 'anchored', 0.5_dp)
  call assert_status(fit, epo_invalid_input)

  singular = reshape([0.04_dp, 0.04_dp, 0.04_dp, 0.04_dp], [2,2])
  fit = epo_from_covariance(singular, signal, 2.0_dp, 'simple', 0.0_dp)
  call assert_status(fit, epo_singular_matrix)

  fit = epo_from_covariance(covariance, [1.0_dp, -2.25_dp], 1.0_dp, &
    'simple', 1.0_dp, normalize=.true.)
  call assert_status(fit, epo_normalization_failure)

  fit = epo_from_covariance(covariance, signal, 2.0_dp, 'anchored', &
    1.0_dp, anchor=anchor, normalize=.false., endogenous=.false.)
  if (.not. fit%ok) error stop 1

  print '(a)', 'test_validation: PASS'

contains

  subroutine assert_status(result, expected)
    type(epo_result), intent(in) :: result
    integer, intent(in) :: expected

    if (result%ok) error stop 1
    if (result%status /= expected) then
      print '(a,2i0)', 'status mismatch: ', result%status, expected
      error stop 1
    end if
  end subroutine assert_status

end program test_validation
