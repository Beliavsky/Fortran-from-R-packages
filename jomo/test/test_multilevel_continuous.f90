program test_multilevel_continuous
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_multilevel, only : jomo1ran_result, jomo1rancon_mcmc
   implicit none

   integer, parameter :: g_count = 6
   integer, parameter :: per_group = 5
   integer, parameter :: n = g_count * per_group
   real(dp) :: y(n, 2)
   logical :: observed(n, 2)
   real(dp) :: x(n, 2)
   real(dp) :: z(n, 1)
   integer :: cluster(n)
   real(dp) :: prior(2, 2)
   real(dp) :: prior_u(2, 2)
   real(dp) :: ug(g_count, 2)
   type(rng_state) :: rng
   type(jomo1ran_result) :: fit
   integer :: i
   integer :: g

   do g = 1, g_count
      ug(g, 1) = 0.20_dp * real(g - 3, dp)
      ug(g, 2) = -0.12_dp * real(g - 3, dp)
   end do
   do i = 1, n
      g = (i - 1) / per_group + 1
      cluster(i) = g
      x(i, 1) = 1.0_dp
      x(i, 2) = real(mod(i - 1, per_group), dp) / real(per_group - 1, dp)
      z(i, 1) = 1.0_dp
      y(i, 1) = 1.0_dp + 0.7_dp * x(i, 2) + ug(g, 1) + 0.05_dp * sin(real(i, dp))
      y(i, 2) = -0.4_dp + 0.3_dp * x(i, 2) + ug(g, 2) + 0.04_dp * cos(real(i, dp))
   end do
   observed = .true.
   observed(4, 1) = .false.
   observed(14, 2) = .false.
   observed(26, :) = .false.
   prior = 0.0_dp
   prior_u = 0.0_dp
   prior(1, 1) = 1.0_dp
   prior(2, 2) = 1.0_dp
   prior_u(1, 1) = 1.0_dp
   prior_u(2, 2) = 1.0_dp

   call rng_seed(rng, 11223344_i8)
   call jomo1rancon_mcmc(rng, y, observed, x, z, cluster, 24, prior, prior_u, fit)

   if (.not. all(ieee_is_finite(fit%level1%continuous))) error stop "non-finite multilevel imputation"
   if (.not. all(ieee_is_finite(fit%random_effects))) error stop "non-finite random effects"
   if (.not. is_spd(fit%level1%omega)) error stop "residual covariance is not SPD"
   if (.not. is_spd(fit%random_covariance)) error stop "random covariance is not SPD"
   if (maxval(abs(pack(fit%level1%continuous - y, observed))) > 1.0e-12_dp) error stop "observed values changed"
   if (size(fit%random_effects, 1) /= g_count) error stop "wrong number of clusters"

   print '(a)', "test_multilevel_continuous: PASS"
end program test_multilevel_continuous
