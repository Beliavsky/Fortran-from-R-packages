! SPDX-License-Identifier: MIT
program simple_epo_example
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

  fit = epo_from_covariance(covariance, signal, 3.0_dp, 'simple', 0.5_dp)
  if (.not. fit%ok) error stop trim(fit%message)

  print '(*(f12.7,1x))', fit%weights
end program simple_epo_example
