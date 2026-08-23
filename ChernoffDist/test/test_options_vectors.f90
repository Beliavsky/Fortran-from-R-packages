program test_options_vectors
  use chernoffdist, only: dp, dchern, pchern, qchern, dchern_vector, pchern_vector, qchern_vector
  implicit none

  real(dp), parameter :: x(3) = [-0.5_dp, 0.0_dp, 0.5_dp]
  real(dp), parameter :: p(3) = [0.1_dp, 0.5_dp, 0.9_dp]
  real(dp) :: dv(3), pv(3), qv(3)
  integer :: i, failures

  failures = 0
  call dchern_vector(x, dv)
  call pchern_vector(x, pv)
  call qchern_vector(p, qv)

  do i = 1, 3
    if (abs(dv(i) - dchern(x(i))) > 1.0e-12_dp) failures = failures + 1
    if (abs(pv(i) - pchern(x(i))) > 1.0e-12_dp) failures = failures + 1
    if (abs(qv(i) - qchern(p(i))) > 1.0e-10_dp) failures = failures + 1
  end do

  if (abs(pchern(0.5_dp, lower_tail=.false.) - (1.0_dp - pchern(0.5_dp))) > 1.0e-12_dp) failures = failures + 1
  if (abs(exp(pchern(0.5_dp, log_p=.true.)) - pchern(0.5_dp)) > 1.0e-12_dp) failures = failures + 1
  if (abs(qchern(log(0.9_dp), log_p=.true.) - qchern(0.9_dp)) > 1.0e-9_dp) failures = failures + 1
  if (abs(qchern(0.1_dp, lower_tail=.false.) - qchern(0.9_dp)) > 1.0e-9_dp) failures = failures + 1

  if (failures /= 0) then
    print '(a,i0)', 'test_options_vectors: FAIL ', failures
    error stop 1
  end if
  print '(a)', 'test_options_vectors: PASS'
end program test_options_vectors
