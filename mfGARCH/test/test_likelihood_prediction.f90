program test_likelihood_prediction
  use mfgarch
  implicit none
  type(mfgarch_model) :: model
  real(dp), allocatable :: contributions(:), tau(:), g(:), residuals(:), forecasts(:)
  real(dp) :: llh, expected, omega
  real(dp) :: returns(5)
  integer :: period(5), horizons(3), status, i

  model%asymmetric = .false.
  model%gamma = 0.0_dp
  model%k = 0
  model%mu = 0.01_dp
  model%alpha = 0.08_dp
  model%beta = 0.87_dp
  model%m = log(1.2_dp)
  returns = [0.10_dp,-0.20_dp,0.05_dp,0.12_dp,-0.08_dp]
  period = [1,2,3,4,5]

  call likelihood_contributions(model, returns, period, contributions, tau, g, residuals, status)
  call assert_true(status == mfgarch_success, 'likelihood status')
  llh = log_likelihood(model, returns, period, status)
  expected = sum(contributions)
  call assert_close(llh, expected, 1.0e-12_dp, 'likelihood sum')
  do i = 1, size(returns)
    call assert_close(residuals(i), (returns(i)-model%mu)/sqrt(tau(i)*g(i)), &
      1.0e-12_dp, 'residual identity')
  end do

  horizons = [1,2,5]
  call predict_variance(model, horizons, exp(model%m), returns(5), g(5), tau(5), &
    forecasts, status)
  call assert_true(status == mfgarch_success, 'prediction status')
  omega = 1.0_dp - model%alpha - model%beta
  expected = exp(model%m) * forecast_garch(omega, model%alpha, model%beta, 0.0_dp, &
    g(5), (returns(5)-model%mu)/sqrt(tau(5)), 1)
  call assert_close(forecasts(1), expected, 1.0e-12_dp, 'one-step prediction')
  call assert_true(all(forecasts > 0.0_dp), 'positive predictions')

  print '(a)', 'test_likelihood_prediction: PASS'

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected_value, tolerance, label)
    real(dp), intent(in) :: actual, expected_value, tolerance
    character(len=*), intent(in) :: label
    call assert_true(abs(actual-expected_value) <= tolerance*(1.0_dp+abs(expected_value)), label)
  end subroutine assert_close

end program test_likelihood_prediction
