program test_tukey_fixture
  use multcomp_kinds, only : dp
  use multcomp_types, only : parm_type, glht_type, mtest_type, confidence_interval_type, contrast_matrix_type
  use multcomp_parm, only : make_parm
  use multcomp_contrasts, only : contr_mat
  use multcomp_glht, only : glht_fit, glht_test, glht_confint
  implicit none

  type(contrast_matrix_type) :: tukey
  type(confidence_interval_type) :: intervals
  type(glht_type) :: hypotheses
  type(mtest_type) :: tests
  type(parm_type) :: parameters
  real(dp) :: beta(3)
  real(dp) :: sigma(3, 3)
  real(dp) :: mean_variance

  beta = [0.0_dp, -10.0_dp, -14.722222_dp]
  mean_variance = 3.872_dp**2 / 2.0_dp
  sigma = 0.0_dp
  sigma(1, 1) = mean_variance
  sigma(2, 2) = mean_variance
  sigma(3, 3) = mean_variance

  call make_parm(beta, sigma, parameters, df=50.0_dp)
  call contr_mat([18.0_dp, 18.0_dp, 18.0_dp], 'Tukey', tukey)
  call glht_fit(parameters, tukey%value, hypotheses)
  call glht_test(hypotheses, 'single-step', tests)
  call glht_confint(hypotheses, 0.95_dp, intervals)
  if (.not. tests%ok .or. .not. intervals%ok) error stop 'warpbreaks-style Tukey fixture failed'

  ! Reference values are from upstream tests/regtest-Tukey.Rout.save.
  if (abs(tests%pvalue(1) - 0.03369_dp) > 0.003_dp) error stop 'Tukey first adjusted p-value mismatch'
  if (abs(tests%pvalue(2) - 0.00105_dp) > 8.0e-4_dp) error stop 'Tukey second adjusted p-value mismatch'
  if (abs(tests%pvalue(3) - 0.44740_dp) > 0.006_dp) error stop 'Tukey third adjusted p-value mismatch'
  if (abs(intervals%critical - 2.4156_dp) > 0.025_dp) error stop 'Tukey simultaneous critical value mismatch'
end program test_tukey_fixture
