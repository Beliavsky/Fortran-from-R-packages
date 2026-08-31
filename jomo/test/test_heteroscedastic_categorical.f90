program test_heteroscedastic_categorical
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_heteroscedastic, only : jomo1ranhr_result, jomo1ranhr_mixed_mcmc
   implicit none

   integer, parameter :: g_count = 5
   integer, parameter :: per_group = 7
   integer, parameter :: n = g_count * per_group
   real(dp) :: y_con(n, 1)
   logical :: con_observed(n, 1)
   integer :: y_cat(n, 1)
   logical :: cat_observed(n, 1)
   integer :: n_levels(1)
   real(dp) :: x(n, 2)
   real(dp) :: z(n, 1)
   integer :: cluster(n)
   real(dp) :: prior(2, 2)
   real(dp) :: prior_u(2, 2)
   type(rng_state) :: rng
   type(jomo1ranhr_result) :: fit
   integer :: i
   integer :: g

   n_levels = 2
   do i = 1, n
      g = (i - 1) / per_group + 1
      cluster(i) = g
      x(i, 1) = 1.0_dp
      x(i, 2) = real(mod(i - 1, per_group), dp) / real(per_group - 1, dp)
      z(i, 1) = 1.0_dp
      y_con(i, 1) = 0.2_dp + 0.5_dp * x(i, 2) + 0.10_dp * real(g - 3, dp) + &
         (0.02_dp + 0.01_dp * real(g, dp)) * cos(real(i, dp))
      if (x(i, 2) + 0.12_dp * real(g - 3, dp) > 0.55_dp) then
         y_cat(i, 1) = 1
      else
         y_cat(i, 1) = 2
      end if
   end do
   con_observed = .true.
   cat_observed = .true.
   con_observed(4, 1) = .false.
   cat_observed(15, 1) = .false.
   con_observed(29, 1) = .false.
   cat_observed(29, 1) = .false.
   prior = 0.0_dp
   prior_u = 0.0_dp
   prior(1, 1) = 1.0_dp
   prior(2, 2) = 1.0_dp
   prior_u(1, 1) = 1.0_dp
   prior_u(2, 2) = 1.0_dp

   call rng_seed(rng, 27182818_i8)
   call jomo1ranhr_mixed_mcmc(rng, y_con, con_observed, y_cat, cat_observed, n_levels, x, z, cluster, 18, &
      prior, prior_u, 5.0_dp, fit, hierarchy_df_prior=2.0_dp)

   if (.not. all(ieee_is_finite(fit%level1%latent))) error stop "heteroscedastic categorical latent state is not finite"
   if (.not. is_spd(fit%random_covariance)) error stop "heteroscedastic categorical random covariance is not SPD"
   if (.not. is_spd(fit%hierarchy_scale)) error stop "heteroscedastic categorical hierarchy scale is not SPD"
   if (.not. ieee_is_finite(fit%hierarchy_df) .or. fit%hierarchy_df < 2.0_dp) &
      error stop "heteroscedastic categorical sampled degrees are invalid"
   do g = 1, g_count
      if (.not. is_spd(fit%cluster_covariance(g, :, :))) error stop "heteroscedastic categorical covariance is not SPD"
      if (abs(fit%cluster_covariance(g, 2, 2) - 1.0_dp) > 1.0e-12_dp) &
         error stop "heteroscedastic binary latent variance changed"
   end do
   do i = 1, n
      if (cat_observed(i, 1)) then
         if (fit%level1%categorical(i, 1) /= y_cat(i, 1)) error stop "observed heteroscedastic category changed"
      end if
   end do

   print '(a)', "test_heteroscedastic_categorical: PASS"
end program test_heteroscedastic_categorical
