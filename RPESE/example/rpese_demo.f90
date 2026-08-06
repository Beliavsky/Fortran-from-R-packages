program rpese_demo
  use rpese, only : dp, rpese_options, se_result, se_if_iid, se_if_cor_pw, &
    mean_se, es_se, rpese_success
  implicit none
  integer, parameter :: n = 96
  real(dp) :: returns(n)
  type(rpese_options) :: options
  type(se_result) :: iid_result, correlated_result, es_result
  integer :: i

  returns(1) = 0.0_dp
  do i = 2, n
    returns(i) = 0.002_dp + 0.45_dp * returns(i - 1) + &
      0.018_dp * sin(0.61_dp * real(i, dp)) + 0.006_dp * cos(1.37_dp * real(i, dp))
  end do

  options = rpese_options()
  options%polynomial_degree = 4
  options%num_lambda = 12
  options%cv_folds = 4
  options%cv_repeats = 1
  options%max_iterations = 400

  call mean_se(returns, iid_result, se_if_iid, options)
  call mean_se(returns, correlated_result, se_if_cor_pw, options)
  call es_se(returns, es_result, confidence=0.95_dp, method=se_if_iid, options=options)

  if (iid_result%status /= rpese_success .or. correlated_result%status /= rpese_success .or. &
      es_result%status /= rpese_success) error stop 'RPESE example failed'

  print '(a,f12.7,a,f12.7)', 'Mean estimate: ', iid_result%estimate, &
    '  iid SE: ', iid_result%standard_error
  print '(a,f12.7,a,f9.5)', 'Prewhitened correlated SE: ', correlated_result%standard_error, &
    '  AR(1): ', correlated_result%ar1_coefficient
  print '(a,f12.7,a,f12.7)', '95% ES estimate: ', es_result%estimate, &
    '  iid SE: ', es_result%standard_error
end program rpese_demo
