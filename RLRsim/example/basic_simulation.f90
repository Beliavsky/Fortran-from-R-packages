program basic_simulation
  use rlrsim_api, only : dp, rlrsim_result, rlrt_sim
  implicit none

  real(dp) :: x(8, 2)
  real(dp) :: z(8, 2)
  real(dp) :: sqrt_sigma(2, 2)
  type(rlrsim_result) :: sim
  integer :: i

  do i = 1, 8
    x(i, 1) = 1.0_dp
    x(i, 2) = real(i, dp) - 4.5_dp
  end do
  z = 0.0_dp
  z(1:4, 1) = 1.0_dp
  z(5:8, 2) = 1.0_dp
  sqrt_sigma = 0.0_dp
  sqrt_sigma(1, 1) = 1.0_dp
  sqrt_sigma(2, 2) = 1.0_dp

  call rlrt_sim(x, z, sqrt_sigma, sim, seed = 2026, nsim = 1000, gridlength = 100)
  print '(a,f10.5)', "Mean simulated RLRT: ", sum(sim%statistic) / real(size(sim%statistic), dp)
  print '(a,f10.5)', "Fraction at zero:     ", &
    real(count(sim%statistic == 0.0_dp), dp) / real(size(sim%statistic), dp)
end program basic_simulation
