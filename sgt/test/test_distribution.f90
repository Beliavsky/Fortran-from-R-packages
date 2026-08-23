program test_distribution
  use sgt
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
  implicit none
  real(dp) :: inf
  integer :: failures
  failures = 0
  inf = ieee_value(0.0_dp, ieee_positive_inf)
  call check(dsgt(0.8_dp, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp), &
    0.28170033592876553_dp, 2.0e-12_dp, failures)
  call check(psgt(0.8_dp, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp), &
    0.6813209779753936_dp, 2.0e-12_dp, failures)
  call check(qsgt(0.37_dp, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp), &
    -0.08149591003676554_dp, 5.0e-12_dp, failures)
  call check(dsgt(-0.2_dp, -0.1_dp, 0.9_dp, -0.3_dp, 1.3_dp, inf), &
    0.4405338895781508_dp, 2.0e-12_dp, failures)
  call check(psgt(-0.2_dp, -0.1_dp, 0.9_dp, -0.3_dp, 1.3_dp, inf), &
    0.3876426013493751_dp, 2.0e-12_dp, failures)
  call check(qsgt(0.37_dp, -0.1_dp, 0.9_dp, -0.3_dp, 1.3_dp, inf), &
    -0.24077704515864484_dp, 5.0e-12_dp, failures)
  call check(dsgt(0.3_dp, 0.1_dp, 0.7_dp, 0.2_dp, inf, 3.0_dp), &
    1.2371791482634837_dp, 2.0e-12_dp, failures)
  call check(psgt(0.3_dp, 0.1_dp, 0.7_dp, 0.2_dp, inf, 3.0_dp), &
    0.7474358296526967_dp, 2.0e-12_dp, failures)
  if (failures /= 0) error stop 1
  print '(a)', 'test_distribution: PASS'
contains
  subroutine check(actual, expected, tol, failures)
    real(dp), intent(in) :: actual, expected, tol
    integer, intent(inout) :: failures
    if (abs(actual - expected) > tol) then
      print *, 'mismatch: ', actual, expected, abs(actual - expected)
      failures = failures + 1
    end if
  end subroutine check
end program test_distribution
