program simultaneous_inference
  use multcomp, only : dp, parm_type, contrast_matrix_type, glht_type, mtest_type, &
    confidence_interval_type, make_parm, contr_mat, glht_fit, glht_test, glht_confint
  implicit none

  type(parm_type) :: parameters
  type(contrast_matrix_type) :: tukey
  type(glht_type) :: hypotheses
  type(mtest_type) :: tests
  type(confidence_interval_type) :: intervals
  real(dp) :: beta(3)
  real(dp) :: covariance(3, 3)
  integer :: i

  beta = [10.0_dp, 12.0_dp, 15.0_dp]
  covariance = 0.0_dp
  covariance(1, 1) = 1.0_dp
  covariance(2, 2) = 1.0_dp
  covariance(3, 3) = 1.0_dp

  call make_parm(beta, covariance, parameters, df=30.0_dp)
  call contr_mat([12.0_dp, 12.0_dp, 12.0_dp], 'Tukey', tukey)
  call glht_fit(parameters, tukey%value, hypotheses)
  call glht_test(hypotheses, 'single-step', tests)
  call glht_confint(hypotheses, 0.95_dp, intervals)

  write (*, '(a,f10.5)') 'simultaneous critical value: ', intervals%critical
  write (*, '(a)') ' estimate       se          t      adjusted-p'
  do i = 1, size(tests%estimate)
    write (*, '(4f12.5)') tests%estimate(i), tests%standard_error(i), &
      tests%statistic(i), tests%pvalue(i)
  end do
end program simultaneous_inference
