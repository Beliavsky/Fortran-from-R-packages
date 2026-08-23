program test_quantile
  use chernoffdist, only: dp, pchern, qchern
  implicit none

  real(dp), parameter :: probs(7) = [0.01_dp, 0.05_dp, 0.10_dp, 0.50_dp, 0.90_dp, 0.95_dp, 0.99_dp]
  real(dp), parameter :: qref(7) = [ &
    -1.1715343421259545_dp, -0.8450811888016130_dp, -0.6642351965043981_dp, &
    0.0_dp, 0.6642351965043982_dp, 0.8450811888016132_dp, 1.1715343421259536_dp ]
  real(dp) :: q
  integer :: i, failures

  failures = 0
  do i = 1, size(probs)
    q = qchern(probs(i))
    if (abs(q - qref(i)) > 3.0e-8_dp) failures = failures + 1
    if (abs(pchern(q) - probs(i)) > 4.0e-9_dp) failures = failures + 1
  end do

  if (failures /= 0) then
    print '(a,i0)', 'test_quantile: FAIL ', failures
    error stop 1
  end if
  print '(a)', 'test_quantile: PASS'
end program test_quantile
