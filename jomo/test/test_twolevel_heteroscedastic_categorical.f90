program test_twolevel_heteroscedastic_categorical
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_twolevel_heteroscedastic, only : jomo2hr_result, jomo2hr_mixed_mcmc
   implicit none

   integer, parameter :: g_count = 6
   integer, parameter :: per_group = 5
   integer, parameter :: n = g_count * per_group
   real(dp) :: y1_con(n, 1)
   logical :: con1_observed(n, 1)
   integer :: y1_cat(n, 1)
   logical :: cat1_observed(n, 1)
   integer :: levels1(1)
   real(dp) :: x1(n, 2)
   real(dp) :: z(n, 1)
   integer :: cluster(n)
   real(dp) :: y2_con(g_count, 0)
   logical :: con2_observed(g_count, 0)
   integer :: y2_cat(g_count, 1)
   logical :: cat2_observed(g_count, 1)
   integer :: levels2(1)
   real(dp) :: x2(g_count, 1)
   real(dp) :: prior1(2, 2)
   real(dp) :: prior_joint(3, 3)
   type(rng_state) :: rng
   type(jomo2hr_result) :: fit
   integer :: i
   integer :: g

   levels1 = 2
   levels2 = 2
   do g = 1, g_count
      x2(g, 1) = 1.0_dp
      y2_cat(g, 1) = merge(1, 2, mod(g, 2) == 0)
   end do
   do i = 1, n
      g = (i - 1) / per_group + 1
      cluster(i) = g
      x1(i, 1) = 1.0_dp
      x1(i, 2) = real(mod(i - 1, per_group), dp) / real(per_group - 1, dp)
      z(i, 1) = 1.0_dp
      y1_con(i, 1) = 0.3_dp + 0.4_dp * x1(i, 2) + merge(0.12_dp, -0.10_dp, y2_cat(g, 1) == 1) + &
         0.03_dp * sin(real(i, dp))
      if (x1(i, 2) + merge(0.15_dp, -0.10_dp, y2_cat(g, 1) == 1) > 0.5_dp) then
         y1_cat(i, 1) = 1
      else
         y1_cat(i, 1) = 2
      end if
   end do
   con1_observed = .true.
   cat1_observed = .true.
   cat2_observed = .true.
   con1_observed(7, 1) = .false.
   cat1_observed(14, 1) = .false.
   con1_observed(26, 1) = .false.
   cat1_observed(26, 1) = .false.
   cat2_observed(5, 1) = .false.
   prior1 = 0.0_dp
   prior1(1, 1) = 1.0_dp
   prior1(2, 2) = 1.0_dp
   prior_joint = 0.0_dp
   prior_joint(1, 1) = 1.0_dp
   prior_joint(2, 2) = 1.0_dp
   prior_joint(3, 3) = 1.0_dp

   call rng_seed(rng, 14142135_i8)
   call jomo2hr_mixed_mcmc(rng, y1_con, con1_observed, y1_cat, cat1_observed, levels1, x1, z, cluster, &
      y2_con, con2_observed, y2_cat, cat2_observed, levels2, x2, 14, prior1, prior_joint, 5.0_dp, fit, &
      hierarchy_df_prior=2.0_dp)

   if (.not. all(ieee_is_finite(fit%base%level1%latent))) error stop "two-level HR categorical latent state is not finite"
   if (.not. is_spd(fit%base%joint_covariance)) error stop "two-level HR categorical joint covariance is not SPD"
   if (abs(fit%base%joint_covariance(3, 3) - 1.0_dp) > 1.0e-12_dp) &
      error stop "level-2 binary latent variance changed"
   do g = 1, g_count
      if (.not. is_spd(fit%cluster_covariance(g, :, :))) error stop "two-level HR categorical cluster covariance is not SPD"
      if (abs(fit%cluster_covariance(g, 2, 2) - 1.0_dp) > 1.0e-12_dp) &
         error stop "level-1 binary latent variance changed"
   end do
   do i = 1, n
      if (cat1_observed(i, 1)) then
         if (fit%base%level1%categorical(i, 1) /= y1_cat(i, 1)) error stop "observed level-1 category changed"
      end if
   end do
   do g = 1, g_count
      if (cat2_observed(g, 1)) then
         if (fit%base%level2%categorical(g, 1) /= y2_cat(g, 1)) error stop "observed level-2 category changed"
      end if
   end do

   print '(a)', "test_twolevel_heteroscedastic_categorical: PASS"
end program test_twolevel_heteroscedastic_categorical
