program test_two_covariates
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use mfgarch
  implicit none
  type(mfgarch_model) :: model
  real(dp), allocatable :: tau(:), weights(:), weights_two(:)
  integer :: period(5), period_two(5), status
  real(dp) :: covariate(3), covariate_two(3), expected_log_tau

  model%asymmetric = .true.
  model%k = 2
  model%k_two = 1
  model%has_second = .true.
  model%m = 0.0_dp
  model%theta = 0.4_dp
  model%w1 = 1.0_dp
  model%w2 = 2.0_dp
  model%theta_two = -0.2_dp
  model%w1_two = 1.0_dp
  model%w2_two = 1.0_dp
  period = [1,2,3,3,3]
  period_two = [1,2,3,3,3]
  covariate = [1.0_dp,2.0_dp,4.0_dp]
  covariate_two = [3.0_dp,5.0_dp,7.0_dp]

  call build_tau(model, period, covariate, tau, status, period_two, covariate_two)
  call assert_true(status == mfgarch_success, 'two-covariate status')
  call assert_true(ieee_is_nan(tau(1)) .and. ieee_is_nan(tau(2)), 'two-covariate leading NaN')
  expected_log_tau = 0.4_dp*(2.0_dp/3.0_dp*2.0_dp + 1.0_dp/3.0_dp*1.0_dp) - 0.2_dp*5.0_dp
  call assert_close(log(tau(3)), expected_log_tau, 1.0e-12_dp, 'two-covariate tau')
  call beta_weights(model%k, model%w1, model%w2, weights, status)
  call beta_weights(model%k_two, model%w1_two, model%w2_two, weights_two, status)
  call assert_close(sum(weights), 1.0_dp, 1.0e-12_dp, 'first weights sum')
  call assert_close(sum(weights_two), 1.0_dp, 1.0e-12_dp, 'second weights sum')
  call assert_true(model_parameter_count(model) == 8, 'parameter count')

  print '(a)', 'test_two_covariates: PASS'

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call assert_true(abs(actual-expected) <= tolerance*(1.0_dp+abs(expected)), label)
  end subroutine assert_close

end program test_two_covariates
