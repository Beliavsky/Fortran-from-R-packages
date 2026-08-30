program test_glht_core
  use multcomp_kinds, only : dp
  use multcomp_types, only : parm_type, glht_type, confidence_interval_type, mtest_type
  use multcomp_parm, only : make_parm
  use multcomp_glht, only : glht_fit, glht_confint, glht_test
  implicit none

  type(parm_type) :: parameters
  type(glht_type) :: hypotheses
  type(confidence_interval_type) :: intervals
  type(mtest_type) :: tests
  real(dp) :: beta(4)
  real(dp) :: sigma(4, 4)
  real(dp) :: k(3, 4)
  real(dp) :: rhs(3)

  beta = [14.8_dp, 12.6667_dp, 7.3333_dp, 13.1333_dp]
  sigma = 0.0_dp
  sigma(1, 1) = 6.7099_dp / 20.0_dp
  sigma(2, 2) = 6.7099_dp / 3.0_dp
  sigma(3, 3) = 6.7099_dp / 3.0_dp
  sigma(4, 4) = 6.7099_dp / 15.0_dp
  k = 0.0_dp
  k(1, 1) = -1.0_dp
  k(1, 2) = 1.0_dp
  k(2, 1) = -1.0_dp
  k(2, 3) = 1.0_dp
  k(3, 1) = -1.0_dp
  k(3, 4) = 1.0_dp
  rhs = 0.0_dp

  call make_parm(beta, sigma, parameters, df=37.0_dp)
  if (.not. parameters%ok) error stop 'make_parm failed'
  call glht_fit(parameters, k, hypotheses, rhs=rhs, alternative='less')
  if (.not. hypotheses%ok) error stop 'glht_fit failed'
  call glht_confint(hypotheses, 0.90_dp, intervals)
  if (.not. intervals%ok) error stop 'glht_confint failed'
  if (abs(intervals%critical - 1.8428_dp) > 0.04_dp) error stop 'critical value regression failed'
  if (abs(intervals%upper(2) - (-4.51131_dp)) > 0.08_dp) error stop 'one-sided upper limit regression failed'

  call glht_test(hypotheses, 'univariate', tests)
  if (.not. tests%ok) error stop 'univariate test failed'
  if (.not. all(tests%pvalue >= 0.0_dp .and. tests%pvalue <= 1.0_dp)) error stop 'invalid p-value'
end program test_glht_core
