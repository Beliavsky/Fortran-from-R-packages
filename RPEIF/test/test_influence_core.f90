program test_influence_core
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rpeif, only : dp, rpeif_options, influence_from_data, rpeif_success, rpeif_numerical_failure
  implicit none
  real(dp) :: returns(10)
  real(dp), allocatable :: values(:), source_values(:), corrected_values(:)
  type(rpeif_options) :: opts
  character(len=16), parameter :: estimators(8) = [character(len=16) :: &
    'mean', 'sd', 'semisd', 'lpm', 'omegaratio', 'sr', 'sor', 'dsr']
  integer :: i, status
  character(len=160) :: message

  returns = [-0.08_dp, -0.03_dp, -0.01_dp, 0.005_dp, 0.01_dp, &
    0.02_dp, 0.03_dp, 0.06_dp, 0.09_dp, 0.12_dp]
  opts%alpha = 0.2_dp
  opts%risk_free = 0.004_dp

  do i = 1, size(estimators)
    call influence_from_data(trim(estimators(i)), returns, returns, values, opts, status, message)
    call assert_true(status == rpeif_success .or. status == rpeif_numerical_failure, trim(estimators(i))//' status')
    call assert_true(all(ieee_is_finite(values)), trim(estimators(i))//' finite')
  end do

  call influence_from_data('mean', returns, returns, values, opts, status)
  call assert_close(maxval(abs(values - (returns - sum(returns) / real(size(returns), dp)))), &
    0.0_dp, 1.0e-14_dp, 'mean IF formula')

  opts%source_compatibility = .true.
  call influence_from_data('omegaratio', returns, returns, source_values, opts, status)
  opts%source_compatibility = .false.
  call influence_from_data('omegaratio', returns, returns, corrected_values, opts, status)
  call assert_true(maxval(abs(source_values - corrected_values)) > 1.0e-6_dp, 'Omega compatibility switch')

  opts%source_compatibility = .true.
  call influence_from_data('sr', returns, returns, source_values, opts, status)
  opts%source_compatibility = .false.
  call influence_from_data('sr', returns, returns, corrected_values, opts, status)
  call assert_true(maxval(abs(source_values - corrected_values)) > 1.0e-8_dp, 'Sharpe compatibility switch')

  print '(a)', 'test_influence_core: PASS'
contains
  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      print '(a,2es24.14)', trim(label)//' failed: ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', trim(label)//' failed'
      error stop 1
    end if
  end subroutine assert_true
end program test_influence_core
