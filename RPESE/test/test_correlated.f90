program test_correlated
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rpese, only : dp, rpese_options, se_result, rpese_success, &
    se_if_cor, se_if_cor_pw, se_if_cor_adapt, fit_exponential, fit_gamma, estimate_se
  implicit none
  integer, parameter :: n = 72
  real(dp) :: x(n)
  type(rpese_options) :: options
  type(se_result) :: result_cor, result_pw, result_adapt, result_gamma
  integer :: i

  x(1) = 0.01_dp
  do i = 2, n
    x(i) = 0.001_dp + 0.65_dp * x(i - 1) + 0.012_dp * sin(0.73_dp * real(i, dp)) + &
      0.004_dp * cos(1.91_dp * real(i, dp))
  end do

  options = rpese_options()
  options%polynomial_degree = 3
  options%num_lambda = 8
  options%cv_folds = 3
  options%cv_repeats = 1
  options%max_iterations = 300
  options%fitting_method = fit_exponential
  call estimate_se(x, 'mean', se_if_cor, result_cor, options)
  call assert_result(result_cor, 'IFcor exponential')
  call assert_true(allocated(result_cor%coefficients), 'IFcor coefficients')

  call estimate_se(x, 'mean', se_if_cor_pw, result_pw, options)
  call assert_result(result_pw, 'IFcorPW exponential')
  call assert_true(abs(result_pw%ar1_coefficient) < 1.0_dp, 'prewhitening coefficient')

  call estimate_se(x, 'mean', se_if_cor_adapt, result_adapt, options)
  call assert_result(result_adapt, 'IFcorAdapt')
  call assert_true(result_adapt%adaptive_weight >= 0.0_dp .and. &
    result_adapt%adaptive_weight <= 1.0_dp, 'adaptive weight')

  options%fitting_method = fit_gamma
  options%num_lambda = 5
  options%max_iterations = 500
  call estimate_se(x, 'sd', se_if_cor, result_gamma, options)
  call assert_result(result_gamma, 'IFcor gamma')

  print '(a)', 'test_correlated: PASS'
contains
  subroutine assert_result(result, label)
    type(se_result), intent(in) :: result
    character(len=*), intent(in) :: label
    if (result%status /= rpese_success) then
      print '(a,i0,2a)', trim(label) // ' status=', result%status, ' message=', trim(result%message)
      error stop 1
    end if
    call assert_true(ieee_is_finite(result%standard_error), trim(label) // ' finite')
    call assert_true(result%standard_error > 0.0_dp, trim(label) // ' positive')
  end subroutine assert_result

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', trim(label) // ' failed'
      error stop 1
    end if
  end subroutine assert_true
end program test_correlated
