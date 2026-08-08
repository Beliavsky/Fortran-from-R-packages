program test_stats
  use rgenoud, only : dp, sample_moments
  implicit none
  real(dp) :: x(4, 1), mu(1), v(1), s(1), k(1)
  x(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  call sample_moments(x, mu, v, s, k)
  if (abs(mu(1) - 2.5_dp) > 1.0e-12_dp) error stop "mean test failed"
  if (abs(v(1) - 1.25_dp) > 1.0e-12_dp) error stop "variance test failed"
  if (abs(s(1)) > 1.0e-12_dp) error stop "skewness test failed"
end program test_stats
