program test_components
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use mfgarch
  implicit none
  type(mfgarch_model) :: model, recovered
  real(dp), allocatable :: weights(:), log_tau(:), tau(:), g(:), raw(:)
  real(dp), parameter :: tolerance = 1.0e-11_dp
  real(dp) :: standardized(3)
  integer :: status, i
  integer :: period(4)
  real(dp) :: covariate(3)

  call beta_weights(3, 1.0_dp, 2.0_dp, weights, status)
  call assert_true(status == mfgarch_success, 'beta_weights status')
  call assert_close(weights(1), 0.5_dp, tolerance, 'weight 1')
  call assert_close(weights(2), 1.0_dp/3.0_dp, tolerance, 'weight 2')
  call assert_close(weights(3), 1.0_dp/6.0_dp, tolerance, 'weight 3')

  call low_frequency_log_tau([1.0_dp,2.0_dp,3.0_dp,4.0_dp], 3, 0.0_dp, 1.0_dp, &
    1.0_dp, 2.0_dp, log_tau, status)
  call assert_true(status == mfgarch_success, 'low_frequency_log_tau status')
  call assert_true(all([(ieee_is_nan(log_tau(i)), i=1,3)]), 'leading tau values are NaN')
  call assert_close(log_tau(4), 7.0_dp/3.0_dp, tolerance, 'log tau value')

  model%k = 2
  model%asymmetric = .false.
  model%gamma = 0.0_dp
  model%m = 0.1_dp
  model%theta = 0.4_dp
  model%w1 = 1.0_dp
  model%w2 = 2.0_dp
  period = [1,2,3,3]
  covariate = [1.0_dp,2.0_dp,4.0_dp]
  call build_tau(model, period, covariate, tau, status)
  call assert_true(status == mfgarch_success, 'build_tau status')
  call assert_true(ieee_is_nan(tau(1)) .and. ieee_is_nan(tau(2)), 'build_tau leading NaN')
  call assert_close(log(tau(3)), 0.1_dp + 0.4_dp*5.0_dp/3.0_dp, tolerance, 'build_tau value')
  call assert_close(tau(4), tau(3), tolerance, 'period expansion')

  standardized = [1.0_dp,-2.0_dp,0.5_dp]
  call calculate_g(standardized, 0.1_dp, 0.7_dp, 0.2_dp, 1.0_dp, g, status)
  call assert_true(status == mfgarch_success, 'calculate_g status')
  call assert_close(g(1), 1.0_dp, tolerance, 'g1')
  call assert_close(g(2), 0.9_dp, tolerance, 'g2')
  call assert_close(g(3), 1.93_dp, tolerance, 'g3')

  model%asymmetric = .true.
  model%gamma = 0.08_dp
  model%alpha = 0.06_dp
  model%beta = 0.88_dp
  model%unrestricted_weights = .true.
  model%w1 = 1.5_dp
  model%w2 = 4.0_dp
  call model_to_raw(model, raw, status)
  call raw_to_model(raw, model, recovered, status)
  call assert_true(status == mfgarch_success, 'model transform status')
  call assert_close(recovered%alpha, model%alpha, 1.0e-10_dp, 'roundtrip alpha')
  call assert_close(recovered%beta, model%beta, 1.0e-10_dp, 'roundtrip beta')
  call assert_close(recovered%gamma, model%gamma, 1.0e-10_dp, 'roundtrip gamma')
  call assert_close(recovered%w1, model%w1, 1.0e-10_dp, 'roundtrip w1')
  call assert_close(recovered%w2, model%w2, 1.0e-10_dp, 'roundtrip w2')

  print '(a)', 'test_components: PASS'

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    call assert_true(abs(actual-expected) <= tol*(1.0_dp+abs(expected)), label)
  end subroutine assert_close

end program test_components
