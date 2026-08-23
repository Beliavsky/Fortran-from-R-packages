program test_options
  use sgt
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_positive_inf
  implicit none
  real(dp) :: inf, f1, f2, p, q
  integer :: failures
  failures = 0
  inf = ieee_value(0.0_dp, ieee_positive_inf)
  f1 = dsgt(0.4_dp, 0.0_dp, 1.2_dp, 0.1_dp, 1.4_dp, 5.0_dp, log_value=.true.)
  f2 = log(dsgt(0.4_dp, 0.0_dp, 1.2_dp, 0.1_dp, 1.4_dp, 5.0_dp))
  call check(f1, f2, 2.0e-14_dp, failures)
  p = psgt(0.7_dp, 0.0_dp, 1.0_dp, -0.2_dp, 1.5_dp, 6.0_dp)
  call check(psgt(0.7_dp, 0.0_dp, 1.0_dp, -0.2_dp, 1.5_dp, 6.0_dp, &
    lower_tail=.false.), 1.0_dp - p, 2.0e-14_dp, failures)
  call check(exp(psgt(0.7_dp, 0.0_dp, 1.0_dp, -0.2_dp, 1.5_dp, 6.0_dp, &
    log_p=.true.)), p, 2.0e-14_dp, failures)
  q = qsgt(log(0.23_dp), 0.0_dp, 1.0_dp, 0.1_dp, 1.8_dp, 7.0_dp, log_p=.true.)
  call check(psgt(q, 0.0_dp, 1.0_dp, 0.1_dp, 1.8_dp, 7.0_dp), 0.23_dp, 2.0e-12_dp, failures)
  f1 = dsgt(0.3_dp, 0.0_dp, 1.1_dp, 0.0_dp, 2.0_dp, inf, &
    var_adj=.false., sigma_multiplier=2.3_dp)
  f2 = dsgt(0.3_dp, 0.0_dp, 1.1_dp * 2.3_dp, 0.0_dp, 2.0_dp, inf, var_adj=.false.)
  call check(f1, f2, 2.0e-14_dp, failures)
  if (.not. ieee_is_nan(dsgt(0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 2.0_dp, inf))) failures = failures + 1
  if (.not. ieee_is_nan(dsgt(0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, inf))) failures = failures + 1
  if (.not. ieee_is_nan(dsgt(0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.2_dp, 4.0_dp))) failures = failures + 1
  if (failures /= 0) error stop 1
  print '(a)', 'test_options: PASS'
contains
  subroutine check(actual, expected, tol, failures)
    real(dp), intent(in) :: actual, expected, tol
    integer, intent(inout) :: failures
    if (abs(actual - expected) > tol) then
      print *, 'mismatch: ', actual, expected
      failures = failures + 1
    end if
  end subroutine check
end program test_options
