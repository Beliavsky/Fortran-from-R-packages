program test_garch_core
  use rumidas
  implicit none
  integer, parameter :: n = 6, k = 2
  real(dp) :: ret(n), mv(k+1,n), param(5), tau, g(n), expected_ll(n), mu, variance
  real(dp), allocatable :: ll(:), cv(:), lr(:), sr(:), wrapped(:)
  type(garch_midas_spec) :: spec
  integer :: i, status

  ret = [0.01_dp, -0.02_dp, 0.015_dp, -0.01_dp, 0.005_dp, -0.012_dp]
  mv = 0.25_dp
  param = [0.10_dp, 0.75_dp, log(0.04_dp), 0.0_dp, 2.0_dp]
  spec%model = RUMIDAS_GM
  spec%skew = .false.
  spec%distribution = RUMIDAS_NORMAL
  spec%lag_function = RUMIDAS_BETA_LAG
  spec%k1 = k

  call garch_midas_evaluate(param, spec, ret, mv, ll, cv, lr, sr, status)
  call check(status == RUMIDAS_SUCCESS, 'core status')
  tau = 0.04_dp
  g(1) = 1.0_dp
  do i = 2, n
    g(i) = (1.0_dp - param(1) - param(2)) + param(1) * ret(i-1)**2 / tau + param(2) * g(i-1)
  end do
  mu = sum(ret) / real(n, dp)
  do i = 1, n
    variance = g(i) * tau
    expected_ll(i) = -0.5_dp * (log(2.0_dp * acos(-1.0_dp)) + log(variance) + (ret(i)-mu)**2/variance)
  end do
  call check(maxval(abs(sr-g)) < 1.0e-13_dp, 'short recursion')
  call check(maxval(abs(lr-sqrt(tau))) < 1.0e-13_dp, 'long run')
  call check(maxval(abs(cv-sqrt(g*tau))) < 1.0e-13_dp, 'conditional')
  call check(maxval(abs(ll-expected_ll)) < 1.0e-12_dp, 'normal likelihood')

  call gm_loglik_no_skew(param, ret, mv, k, wrapped, status)
  call check(status == RUMIDAS_SUCCESS, 'wrapper status')
  call check(maxval(abs(wrapped-ll)) < 1.0e-13_dp, 'wrapper identity')

  print '(a)', 'test_garch_core: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) error stop message
  end subroutine check
end program test_garch_core
