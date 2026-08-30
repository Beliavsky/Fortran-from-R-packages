program test_adjustments
  use multcomp_kinds, only : dp
  use multcomp_types, only : parm_type, glht_type, mtest_type, global_test_type, cld_type
  use multcomp_math, only : p_adjust
  use multcomp_parm, only : make_parm
  use multcomp_contrasts, only : contr_mat
  use multcomp_types, only : contrast_matrix_type
  use multcomp_glht, only : glht_fit, glht_test, glht_global_test, glht_coefficients
  use multcomp_helpers, only : cld_from_tukey_test
  implicit none

  type(contrast_matrix_type) :: tukey
  type(global_test_type) :: global
  type(cld_type) :: letters
  type(glht_type) :: coefficient_tests
  type(glht_type) :: hypotheses
  type(mtest_type) :: adjusted_test
  type(mtest_type) :: shaffer
  type(mtest_type) :: univariate_test
  type(mtest_type) :: westfall
  type(parm_type) :: parameters
  real(dp), allocatable :: adjusted(:)
  real(dp) :: beta(3)
  real(dp) :: expected(5)
  real(dp) :: p(5)
  real(dp) :: sigma(3, 3)
  logical :: ok

  p = [0.01_dp, 0.04_dp, 0.03_dp, 0.002_dp, 0.5_dp]

  expected = [0.05_dp, 0.20_dp, 0.15_dp, 0.01_dp, 1.0_dp]
  call p_adjust(p, 'bonferroni', adjusted, ok)
  call check_vector(adjusted, expected, 1.0e-14_dp, 'Bonferroni')

  expected = [0.04_dp, 0.09_dp, 0.09_dp, 0.01_dp, 0.5_dp]
  call p_adjust(p, 'holm', adjusted, ok)
  call check_vector(adjusted, expected, 1.0e-14_dp, 'Holm')

  expected = [0.04_dp, 0.08_dp, 0.08_dp, 0.01_dp, 0.5_dp]
  call p_adjust(p, 'hochberg', adjusted, ok)
  call check_vector(adjusted, expected, 1.0e-14_dp, 'Hochberg')

  expected = [0.04_dp, 0.08_dp, 0.06_dp, 0.01_dp, 0.5_dp]
  call p_adjust(p, 'hommel', adjusted, ok)
  call check_vector(adjusted, expected, 1.0e-14_dp, 'Hommel')

  expected = [0.025_dp, 0.05_dp, 0.05_dp, 0.01_dp, 0.5_dp]
  call p_adjust(p, 'BH', adjusted, ok)
  call check_vector(adjusted, expected, 1.0e-14_dp, 'BH')

  expected = [0.057083333333333333_dp, 0.11416666666666667_dp, &
    0.11416666666666667_dp, 0.022833333333333333_dp, 1.0_dp]
  call p_adjust(p, 'BY', adjusted, ok)
  call check_vector(adjusted, expected, 2.0e-14_dp, 'BY')

  beta = [0.0_dp, 0.7_dp, 1.6_dp]
  sigma = 0.0_dp
  sigma(1, 1) = 0.09_dp
  sigma(2, 2) = 0.09_dp
  sigma(3, 3) = 0.09_dp
  call make_parm(beta, sigma, parameters, df=40.0_dp)
  if (.not. parameters%ok) error stop 'parameter construction failed'
  call contr_mat([10.0_dp, 10.0_dp, 10.0_dp], 'Tukey', tukey)
  if (.not. tukey%ok) error stop 'Tukey contrast construction failed'
  call glht_fit(parameters, tukey%value, hypotheses)
  if (.not. hypotheses%ok) error stop 'Tukey GLHT construction failed'

  call glht_test(hypotheses, 'univariate', univariate_test)
  call glht_test(hypotheses, 'single-step', adjusted_test)
  call glht_test(hypotheses, 'Shaffer', shaffer)
  call glht_test(hypotheses, 'Westfall', westfall)
  if (.not. adjusted_test%ok .or. .not. shaffer%ok .or. .not. westfall%ok) &
    error stop 'simultaneous adjustment failed'
  if (any(adjusted_test%pvalue + 1.0e-10_dp < univariate_test%pvalue)) &
    error stop 'single-step p-values must dominate marginal p-values'
  if (any(shaffer%pvalue + 1.0e-10_dp < univariate_test%pvalue)) &
    error stop 'Shaffer p-values must dominate marginal p-values'
  if (any(westfall%pvalue + 1.0e-10_dp < univariate_test%pvalue)) &
    error stop 'Westfall p-values must dominate marginal p-values'
  if (any(westfall%pvalue > shaffer%pvalue + 0.005_dp)) &
    error stop 'Westfall unexpectedly exceeds Shaffer on the three-group fixture'

  call glht_global_test(hypotheses, global, request_f=.true.)
  if (.not. global%ok .or. .not. global%f_test) error stop 'global F test failed'
  if (global%rank /= 2) error stop 'Tukey global hypothesis rank should be two'
  if (global%pvalue < 0.0_dp .or. global%pvalue > 1.0_dp) error stop 'invalid global p-value'

  call glht_coefficients(parameters, [1, 3], coefficient_tests)
  if (.not. coefficient_tests%ok .or. size(coefficient_tests%estimate) /= 2) &
    error stop 'selected-coefficient cftest equivalent failed'
  if (maxval(abs(coefficient_tests%estimate - [0.0_dp, 1.6_dp])) > 1.0e-14_dp) &
    error stop 'selected-coefficient estimates are wrong'

  call cld_from_tukey_test(shaffer, 3, 0.05_dp, letters)
  if (.not. letters%ok .or. size(letters%letter_matrix, 1) /= 3) &
    error stop 'Tukey compact-letter convenience failed'

contains

  subroutine check_vector(actual, target, tolerance, label)
    real(dp), intent(in) :: actual(:) !! Computed p-values to compare with the reference vector.
    real(dp), intent(in) :: target(:) !! Reference adjusted p-values in matching order.
    real(dp), intent(in) :: tolerance !! Maximum accepted absolute componentwise discrepancy.
    character(len=*), intent(in) :: label !! Adjustment name reported when the regression check fails.

    if (.not. all(abs(actual - target) <= tolerance)) then
      write (*, '(a,1x,*(es16.8,1x))') trim(label)//' actual:', actual
      write (*, '(a,1x,*(es16.8,1x))') trim(label)//' target:', target
      error stop 'p-adjust regression mismatch'
    end if
  end subroutine check_vector

end program test_adjustments
