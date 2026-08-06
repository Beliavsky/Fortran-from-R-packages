program test_regression
  use robstattm, only : dp, robstattm_control, regression_result, linear_test_result, &
    model_selection_result, &
    lmrobdet_control, lmrobdet_mm, lmrobdet_dcml, least_squares_fit, &
    lmrobdet_mm_rfpe, lmrobdet_lin_test, stepwise_rfpe
  implicit none
  integer, parameter :: n = 48, p = 3
  real(dp) :: x(n, p), xr(n, 2), y(n), t, rfpe, minimum_rho, penalty
  type(robstattm_control) :: control
  type(regression_result) :: robust_fit, ls_fit, dcml, restricted_fit
  type(linear_test_result) :: test
  type(model_selection_result) :: selection
  integer :: i

  do i = 1, n
    t = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
    x(i, :) = [1.0_dp, t, sin(3.0_dp * t)]
    y(i) = 1.5_dp + 2.0_dp * t - 0.7_dp * sin(3.0_dp * t) + 0.08_dp * cos(11.0_dp * t)
  end do
  y(5) = y(5) + 12.0_dp
  y(17) = y(17) - 10.0_dp
  y(39) = y(39) + 8.0_dp

  control = lmrobdet_control(n_resample=120, max_iter=100, tolerance=1.0e-7_dp)
  call lmrobdet_mm(x, y, robust_fit, control)
  call assert_true(allocated(robust_fit%coefficients), 'MM coefficients allocated')
  call assert_true(size(robust_fit%coefficients) == p, 'MM coefficient count')
  call assert_true(abs(robust_fit%coefficients(2) - 2.0_dp) < 0.6_dp, 'MM slope recovery')
  call assert_true(robust_fit%scale > 0.0_dp, 'MM scale')

  call least_squares_fit(x, y, ls_fit)
  call assert_true(allocated(ls_fit%covariance), 'LS covariance')
  call lmrobdet_dcml(x, y, dcml, control)
  call assert_true(allocated(dcml%coefficients), 'DCML coefficients allocated')
  call assert_true(dcml%mixing >= 0.0_dp .and. dcml%mixing <= 1.0_dp, 'DCML mixing')

  rfpe = lmrobdet_mm_rfpe(robust_fit, control, minimum_rho=minimum_rho, penalty=penalty)
  call assert_true(rfpe < huge(1.0_dp), 'finite RFPE')
  call assert_true(minimum_rho >= 0.0_dp .and. penalty >= 0.0_dp, 'RFPE terms')

  xr = x(:, 1:2)
  call lmrobdet_mm(xr, y, restricted_fit, control)
  call lmrobdet_lin_test(robust_fit, restricted_fit, test, control)
  call assert_true(test%df1 == 1, 'linear test numerator degrees')
  call assert_true(test%chi_square_p_value >= 0.0_dp .and. test%chi_square_p_value <= 1.0_dp, &
    'linear test p-value')
  call stepwise_rfpe(x, y, selection, control, direction='backward', max_steps=2)
  call assert_true(allocated(selection%selected_columns), 'RFPE selected columns')
  call assert_true(size(selection%selected_columns) >= 1, 'RFPE nonempty model')
  call assert_true(selection%criterion < huge(1.0_dp), 'RFPE selection criterion')
  print '(a)', 'test_regression: PASS'
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_regression
