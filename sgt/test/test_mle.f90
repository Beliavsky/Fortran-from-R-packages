program test_mle
  use sgt
  implicit none
  integer, parameter :: n = 160
  real(dp) :: data(n), u
  type(sgt_params) :: start
  type(sgt_mle_result) :: fit
  logical :: free(5)
  integer :: i, failures
  failures = 0
  do i = 1, n
    u = (real(i,dp) - 0.5_dp) / real(n,dp)
    data(i) = qsgt(u, 0.7_dp, 1.2_dp, -0.2_dp, 1.6_dp, 5.0_dp)
  end do
  start%mu = 0.4_dp
  start%sigma = 1.0_dp
  start%lambda = -0.05_dp
  start%p = 1.6_dp
  start%q = 5.0_dp
  free = [.true., .true., .true., .false., .false.]
  call sgt_mle_constant(data, start, fit, free=free, max_iter=500)
  if (fit%convcode /= 0) failures = failures + 1
  if (abs(fit%estimate(1) - 0.7_dp) > 0.015_dp) failures = failures + 1
  if (abs(fit%estimate(2) - 1.2_dp) > 0.015_dp) failures = failures + 1
  if (abs(fit%estimate(3) + 0.2_dp) > 0.015_dp) failures = failures + 1
  if (maxval(abs(fit%gradient)) > 2.0e-4_dp) failures = failures + 1
  if (.not. all(fit%std_error > 0.0_dp)) failures = failures + 1
  if (failures /= 0) then
    print *, fit%estimate
    print *, fit%gradient
    print *, fit%std_error
    error stop 1
  end if
  print '(a)', 'test_mle: PASS'
end program test_mle
