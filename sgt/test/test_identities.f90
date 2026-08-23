program test_identities
  use sgt
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
  implicit none
  real(dp) :: inf, x, h, deriv, expected, pi
  real(dp) :: probs(7), qv(7), back(7)
  integer :: i, failures
  failures = 0
  inf = ieee_value(0.0_dp, ieee_positive_inf)
  pi = acos(-1.0_dp)
  call check(dsgt(1.2_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, inf), &
    exp(-0.5_dp * 1.2_dp**2) / sqrt(2.0_dp * pi), 5.0e-13_dp, failures)
  call check(psgt(1.2_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, inf), &
    0.5_dp * erfc(-1.2_dp / sqrt(2.0_dp)), 5.0e-13_dp, failures)
  call check(qsgt(0.8_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, inf), &
    0.8416212335729143_dp, 5.0e-12_dp, failures)
  expected = 1.0_dp / (pi * 1.3_dp * (1.0_dp + ((0.7_dp - 1.0_dp) / 1.3_dp)**2))
  call check(dsgt(0.7_dp, 1.0_dp, 1.3_dp, 0.0_dp, 2.0_dp, 0.5_dp, &
    mean_cent=.false., var_adj=.false., sigma_multiplier=sqrt(2.0_dp)), expected, &
    2.0e-12_dp, failures)
  expected = exp(-abs(0.2_dp - 1.2_dp) / 1.8_dp) / (2.0_dp * 1.8_dp)
  call check(dsgt(0.2_dp, 1.2_dp, 1.8_dp, 0.0_dp, 1.0_dp, inf, &
    mean_cent=.false., var_adj=.false.), expected, 2.0e-12_dp, failures)
  call check(dsgt(1.7_dp, 1.9_dp, 0.7_dp, 0.0_dp, inf, 2.0_dp, &
    var_adj=.false.), 1.0_dp / 1.4_dp, 2.0e-14_dp, failures)
  call check(psgt(qsgt((1.0_dp - 0.25_dp) / 2.0_dp, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp), &
    0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp), (1.0_dp - 0.25_dp) / 2.0_dp, &
    2.0e-12_dp, failures)
  probs = [1.0e-8_dp, 0.01_dp, 0.1_dp, 0.37_dp, 0.8_dp, 0.99_dp, 1.0_dp - 1.0e-8_dp]
  qv = qsgt(probs, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp)
  back = psgt(qv, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp)
  do i = 1, size(probs)
    call check(back(i), probs(i), 2.0e-12_dp, failures)
  end do
  x = 0.55_dp
  h = 1.0e-5_dp
  deriv = (psgt(x + h, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp) - &
    psgt(x - h, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp)) / (2.0_dp * h)
  call check(deriv, dsgt(x, 0.4_dp, 1.3_dp, 0.25_dp, 1.7_dp, 4.2_dp), &
    2.0e-9_dp, failures)
  if (failures /= 0) error stop 1
  print '(a)', 'test_identities: PASS'
contains
  subroutine check(actual, expected, tol, failures)
    real(dp), intent(in) :: actual, expected, tol
    integer, intent(inout) :: failures
    if (abs(actual - expected) > tol) then
      print *, 'mismatch: ', actual, expected, abs(actual - expected)
      failures = failures + 1
    end if
  end subroutine check
end program test_identities
