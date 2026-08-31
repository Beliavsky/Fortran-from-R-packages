program jomo1_continuous_example
   use jomo, only : dp, i8, rng_state, rng_seed, jomo1_result, jomo1con_mcmc
   implicit none

   integer, parameter :: n = 24
   real(dp) :: y(n, 2)
   logical :: observed(n, 2)
   real(dp) :: x(n, 2)
   real(dp) :: prior(2, 2)
   type(rng_state) :: rng
   type(jomo1_result) :: fit
   integer :: i

   do i = 1, n
      x(i, 1) = 1.0_dp
      x(i, 2) = real(i - 1, dp) / real(n - 1, dp)
      y(i, 1) = 0.5_dp + 0.7_dp * x(i, 2) + 0.08_dp * sin(real(i, dp))
      y(i, 2) = -0.3_dp + 0.4_dp * x(i, 2) + 0.06_dp * cos(real(i, dp))
   end do
   observed = .true.
   observed(5, 1) = .false.
   observed(14, 2) = .false.
   observed(20, :) = .false.
   prior = 0.0_dp
   prior(1, 1) = 1.0_dp
   prior(2, 2) = 1.0_dp

   call rng_seed(rng, 20260830_i8)
   call jomo1con_mcmc(rng, y, observed, x, 50, prior, fit)

   print '(a,2f12.6)', "posterior mean beta, intercepts: ", fit%beta_mean(1, :)
   print '(a,2f12.6)', "imputed row 20:                 ", fit%continuous(20, :)
end program jomo1_continuous_example
