program test_cholesterol_fixture
  use multcomp_kinds, only : dp
  use multcomp_types, only : parm_type, glht_type, mtest_type, contrast_matrix_type
  use multcomp_parm, only : make_parm
  use multcomp_contrasts, only : contr_mat
  use multcomp_glht, only : glht_fit, glht_test
  implicit none

  type(contrast_matrix_type) :: tukey
  type(glht_type) :: hypotheses
  type(mtest_type) :: shaffer
  type(mtest_type) :: univariate
  type(mtest_type) :: westfall
  type(parm_type) :: parameters
  real(dp) :: beta(5)
  real(dp) :: sigma(5, 5)
  real(dp) :: mean_variance

  ! These group means reproduce the rounded estimates printed by the upstream
  ! cholesterol example; equal group variances reproduce its 1.443 standard errors.
  beta = [0.0_dp, 3.443_dp, 6.593_dp, 9.579_dp, 15.166_dp]
  mean_variance = 1.443_dp**2 / 2.0_dp
  sigma = 0.0_dp
  sigma(1, 1) = mean_variance
  sigma(2, 2) = mean_variance
  sigma(3, 3) = mean_variance
  sigma(4, 4) = mean_variance
  sigma(5, 5) = mean_variance

  call make_parm(beta, sigma, parameters, df=45.0_dp)
  call contr_mat([10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp], 'Tukey', tukey)
  call glht_fit(parameters, tukey%value, hypotheses)
  if (.not. hypotheses%ok) error stop 'cholesterol fixture GLHT construction failed'

  call glht_test(hypotheses, 'univariate', univariate)
  call glht_test(hypotheses, 'Shaffer', shaffer)
  call glht_test(hypotheses, 'Westfall', westfall)
  if (.not. univariate%ok .or. .not. shaffer%ok .or. .not. westfall%ok) &
    error stop 'cholesterol fixture testing failed'

  ! Upstream multcomp-Ex.Rout.save prints these values rounded to six significant digits.
  if (abs(univariate%pvalue(1) - 0.021333_dp) > 8.0e-5_dp) &
    error stop 'cholesterol univariate p-value mismatch'
  if (abs(shaffer%pvalue(1) - 0.042666_dp) > 1.5e-4_dp) &
    error stop 'cholesterol Shaffer first p-value mismatch'
  if (abs(shaffer%pvalue(5) - 0.042666_dp) > 1.5e-4_dp) &
    error stop 'cholesterol Shaffer fifth p-value mismatch'
  if (abs(shaffer%pvalue(8) - 0.044316_dp) > 1.5e-4_dp) &
    error stop 'cholesterol Shaffer eighth p-value mismatch'
  if (abs(westfall%pvalue(1) - 0.0420_dp) > 2.0e-3_dp) &
    error stop 'cholesterol Westfall first p-value mismatch'
  if (abs(westfall%pvalue(5) - 0.0420_dp) > 2.0e-3_dp) &
    error stop 'cholesterol Westfall fifth p-value mismatch'
  if (abs(westfall%pvalue(8) - 0.0443_dp) > 2.0e-3_dp) &
    error stop 'cholesterol Westfall eighth p-value mismatch'
end program test_cholesterol_fixture
