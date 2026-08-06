program test_bootstrap_api
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rpese, only : dp, rpese_options, se_result, rpese_success, &
    se_boot_iid, se_boot_cor, mean_se, es_se, var_se, lpm_se, omegaratio_se, &
    se_method_from_name, fitting_method_from_name, frequency_mode_from_name, &
    fit_exponential, frequency_decimate
  implicit none
  integer, parameter :: n = 50
  real(dp) :: x(n)
  type(rpese_options) :: options
  type(se_result) :: result
  integer :: i

  do i = 1, n
    x(i) = 0.004_dp + 0.02_dp * sin(0.47_dp * real(i, dp)) + 0.006_dp * cos(1.31_dp * real(i, dp))
  end do
  options = rpese_options()
  options%bootstrap_replicates = 80
  options%seed = 2468
  call mean_se(x, result, se_boot_iid, options)
  call assert_result(result, 'iid bootstrap')
  call mean_se(x, result, se_boot_cor, options)
  call assert_result(result, 'block bootstrap')

  call es_se(x, result, confidence=0.9_dp, method=se_boot_iid, options=options)
  call assert_result(result, 'ES wrapper')
  call var_se(x, result, confidence=0.9_dp, method=se_boot_iid, options=options)
  call assert_result(result, 'VaR wrapper')
  call lpm_se(x, result, threshold=0.0_dp, order=2, method=se_boot_iid, options=options)
  call assert_result(result, 'LPM wrapper')
  call omegaratio_se(x, result, threshold=0.0_dp, method=se_boot_iid, options=options)
  call assert_result(result, 'Omega wrapper')

  call assert_true(se_method_from_name('IFcorPW') > 0, 'method parser')
  call assert_true(fitting_method_from_name('Exponential') == fit_exponential, 'fit parser')
  call assert_true(frequency_mode_from_name('Decimate') == frequency_decimate, 'frequency parser')

  print '(a)', 'test_bootstrap_api: PASS'
contains
  subroutine assert_result(result, label)
    type(se_result), intent(in) :: result
    character(len=*), intent(in) :: label
    if (result%status /= rpese_success) then
      print '(a,i0,2a)', trim(label) // ' status=', result%status, ' message=', trim(result%message)
      error stop 1
    end if
    call assert_true(ieee_is_finite(result%standard_error), trim(label) // ' finite')
    call assert_true(result%standard_error >= 0.0_dp, trim(label) // ' nonnegative')
  end subroutine assert_result

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', trim(label) // ' failed'
      error stop 1
    end if
  end subroutine assert_true
end program test_bootstrap_api
