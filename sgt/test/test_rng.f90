program test_rng
  use sgt
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
  implicit none
  integer, parameter :: n = 50000
  real(dp), allocatable :: x(:)
  real(dp) :: inf, meanx, sdx
  integer, allocatable :: seed(:)
  integer :: nseed, i, failures
  failures = 0
  inf = ieee_value(0.0_dp, ieee_positive_inf)
  call random_seed(size=nseed)
  allocate(seed(nseed), x(n))
  do i = 1, nseed
    seed(i) = 137 + 17 * i
  end do
  call random_seed(put=seed)
  call rsgt(x, mu=0.8_dp, sigma=1.4_dp, lambda=-0.25_dp, p=1.7_dp, q=6.0_dp)
  meanx = sum(x) / real(n, dp)
  sdx = sqrt(sum((x - meanx)**2) / real(n - 1, dp))
  if (abs(meanx - 0.8_dp) > 0.035_dp) failures = failures + 1
  if (abs(sdx - 1.4_dp) > 0.035_dp) failures = failures + 1
  call random_seed(put=seed)
  call rsgt(x, mu=-0.3_dp, sigma=0.9_dp, p=2.0_dp, q=inf)
  meanx = sum(x) / real(n, dp)
  sdx = sqrt(sum((x - meanx)**2) / real(n - 1, dp))
  if (abs(meanx + 0.3_dp) > 0.025_dp) failures = failures + 1
  if (abs(sdx - 0.9_dp) > 0.025_dp) failures = failures + 1
  if (failures /= 0) then
    print *, meanx, sdx
    error stop 1
  end if
  print '(a)', 'test_rng: PASS'
end program test_rng
