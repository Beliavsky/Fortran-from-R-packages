! SPDX-License-Identifier: MIT
program epo_demo
  use epo, only : dp, epo_optimize, epo_result
  implicit none

  real(dp) :: anchor(3), returns(8,3), signal(3)
  type(epo_result) :: fit

  returns = reshape([ &
    0.012_dp, -0.004_dp, 0.009_dp, 0.015_dp, -0.007_dp, 0.006_dp, 0.011_dp, 0.003_dp, &
    0.006_dp,  0.010_dp, 0.002_dp, 0.013_dp, -0.003_dp, 0.008_dp, 0.004_dp, 0.009_dp, &
   -0.002_dp,  0.005_dp, 0.011_dp, 0.004_dp,  0.009_dp, 0.001_dp, 0.007_dp, 0.010_dp  &
  ], [8,3])
  signal = [0.08_dp, 0.06_dp, 0.05_dp]
  anchor = [0.4_dp, 0.3_dp, 0.3_dp]

  fit = epo_optimize(returns, signal, 10.0_dp, 'anchored', 0.5_dp, &
    anchor=anchor)
  if (.not. fit%ok) error stop trim(fit%message)

  print '(a)', 'Anchored EPO weights:'
  print '(*(f12.7,1x))', fit%weights
  print '(a,f12.7)', 'Endogenous gamma: ', fit%gamma
end program epo_demo
