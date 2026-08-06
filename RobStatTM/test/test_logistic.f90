program test_logistic
  use robstattm, only : dp, logistic_result, logreg_by, logreg_wby, logreg_wml
  implicit none
  integer, parameter :: n = 64
  real(dp) :: x(n, 2), y(n), t, probability, sequence
  type(logistic_result) :: by, wby, wml
  integer :: i

  do i = 1, n
    t = -2.5_dp + 5.0_dp * real(i - 1, dp) / real(n - 1, dp)
    x(i, 1) = t
    x(i, 2) = cos(1.7_dp * t)
    probability = 1.0_dp / (1.0_dp + exp(-(-0.25_dp + 0.9_dp * t - 0.55_dp * x(i, 2))))
    sequence = real(mod(37 * i + 11, 101), dp) / 101.0_dp
    y(i) = merge(1.0_dp, 0.0_dp, sequence < probability)
  end do
  x(n, :) = [18.0_dp, -12.0_dp]
  y(n) = 0.0_dp

  call logreg_by(x, y, by, max_iter=300)
  call assert_logistic(by, n, 3, 'BY')
  call logreg_wby(x, y, wby, max_iter=300)
  call assert_logistic(wby, n, 3, 'WBY')
  call logreg_wml(x, y, wml, max_iter=100)
  call assert_logistic(wml, n, 3, 'WML')
  call assert_true(wby%leverage_weights(n) < 0.5_dp .or. wml%leverage_weights(n) < 0.5_dp, &
    'high-leverage point filtered')
  print '(a)', 'test_logistic: PASS'
contains
  subroutine assert_logistic(fit, observations, parameters, name)
    type(logistic_result), intent(in) :: fit
    integer, intent(in) :: observations, parameters
    character(len=*), intent(in) :: name
    call assert_true(allocated(fit%coefficients), trim(name) // ' coefficients')
    call assert_true(size(fit%coefficients) == parameters, trim(name) // ' parameter count')
    call assert_true(size(fit%fitted) == observations, trim(name) // ' fitted count')
    call assert_true(all(fit%fitted > 0.0_dp .and. fit%fitted < 1.0_dp), trim(name) // ' probabilities')
  end subroutine assert_logistic

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_logistic
