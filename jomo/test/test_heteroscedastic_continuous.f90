program test_heteroscedastic_continuous
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_heteroscedastic, only : jomo1ranhr_result, jomo1ranconhr_mcmc
   implicit none

   integer, parameter :: g_count = 5
   integer, parameter :: per_group = 6
   integer, parameter :: n = g_count * per_group
   real(dp) :: y(n, 2)
   logical :: observed(n, 2)
   real(dp) :: x(n, 2)
   real(dp) :: z(n, 1)
   integer :: cluster(n)
   real(dp) :: prior(2, 2)
   real(dp) :: prior_u(2, 2)
   type(rng_state) :: rng
   type(jomo1ranhr_result) :: fit
   integer :: i
   integer :: g

   do i = 1, n
      g = (i - 1) / per_group + 1
      cluster(i) = g
      x(i, 1) = 1.0_dp
      x(i, 2) = real(mod(i - 1, per_group), dp) / real(per_group - 1, dp)
      z(i, 1) = 1.0_dp
      y(i, 1) = 0.5_dp + 0.4_dp * x(i, 2) + 0.10_dp * real(g - 3, dp) + &
         (0.02_dp + 0.01_dp * real(g, dp)) * sin(real(i, dp))
      y(i, 2) = -0.2_dp + 0.3_dp * x(i, 2) - 0.08_dp * real(g - 3, dp) + &
         (0.03_dp + 0.008_dp * real(g, dp)) * cos(real(i, dp))
   end do
   observed = .true.
   observed(3, 1) = .false.
   observed(17, 2) = .false.
   observed(28, :) = .false.
   prior = 0.0_dp
   prior_u = 0.0_dp
   prior(1, 1) = 1.0_dp
   prior(2, 2) = 1.0_dp
   prior_u(1, 1) = 1.0_dp
   prior_u(2, 2) = 1.0_dp

   call rng_seed(rng, 55443322_i8)
   call jomo1ranconhr_mcmc(rng, y, observed, x, z, cluster, 16, prior, prior_u, 4.0_dp, fit)

   if (.not. all(ieee_is_finite(fit%level1%continuous))) error stop "non-finite heterogeneous imputation"
   if (.not. all(ieee_is_finite(fit%random_effects))) error stop "non-finite heterogeneous random effects"
   if (.not. is_spd(fit%random_covariance)) error stop "heterogeneous random covariance is not SPD"
   if (.not. is_spd(fit%hierarchy_scale)) error stop "heterogeneous hierarchy scale is not SPD"
   if (.not. ieee_is_finite(fit%hierarchy_df) .or. fit%hierarchy_df < 2.0_dp) &
      error stop "heterogeneous hierarchy degrees update is invalid"
   do g = 1, g_count
      if (.not. is_spd(fit%cluster_covariance(g, :, :))) error stop "cluster residual covariance is not SPD"
   end do
   if (maxval(abs(pack(fit%level1%continuous - y, observed))) > 1.0e-12_dp) error stop "observed values changed"

   print '(a)', "test_heteroscedastic_continuous: PASS"
end program test_heteroscedastic_continuous
