program test_identities
  use chernoffdist, only: dp, dchern, pchern
  implicit none

  real(dp), parameter :: x(5) = [0.1_dp, 0.3_dp, 0.8_dp, 1.2_dp, 2.0_dp]
  real(dp) :: h, deriv
  integer :: i, failures

  failures = 0
  do i = 1, size(x)
    if (abs(dchern(x(i)) - dchern(-x(i))) > 2.0e-10_dp) failures = failures + 1
    if (abs(pchern(x(i)) + pchern(-x(i)) - 1.0_dp) > 3.0e-9_dp) failures = failures + 1
  end do

  h = 2.0e-5_dp
  deriv = (pchern(0.5_dp + h) - pchern(0.5_dp - h)) / (2.0_dp * h)
  if (abs(deriv - dchern(0.5_dp)) > 2.0e-6_dp) failures = failures + 1

  if (failures /= 0) then
    print '(a,i0)', 'test_identities: FAIL ', failures
    error stop 1
  end if
  print '(a)', 'test_identities: PASS'
end program test_identities
