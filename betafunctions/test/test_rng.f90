program test_rng
  use betafunctions
  implicit none
  real(dp), allocatable :: x(:)
  integer, allocatable :: ix(:)
  integer, allocatable :: seed(:)
  integer :: nseed

  call random_seed(size=nseed)
  allocate(seed(nseed))
  seed = 1729
  call random_seed(put=seed)

  allocate(x(20000))
  call beta4_random(x, 0.25_dp, 0.75_dp, 5.0_dp, 3.0_dp)
  if (abs(sum(x)/real(size(x),dp) - 0.5625_dp) > 0.004_dp) error stop 1
  if (minval(x) < 0.25_dp .or. maxval(x) > 0.75_dp) error stop 1

  allocate(ix(2000))
  call beta_binomial_random(ix, 20, 0.0_dp, 1.0_dp, 2.5_dp, 4.0_dp)
  if (minval(ix) < 0 .or. maxval(ix) > 20) error stop 1

  print '(a)', 'test_rng: PASS'
end program test_rng
