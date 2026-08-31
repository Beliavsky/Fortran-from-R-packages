program test_twolevel_heteroscedastic_continuous
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_twolevel_heteroscedastic, only : jomo2hr_result, jomo2conhr_mcmc
   implicit none

   integer, parameter :: g_count = 6
   integer, parameter :: per_group = 6
   integer, parameter :: n = g_count * per_group
   real(dp) :: y1(n, 2)
   logical :: observed1(n, 2)
   real(dp) :: x1(n, 2)
   real(dp) :: z(n, 1)
   integer :: cluster(n)
   real(dp) :: y2(g_count, 1)
   logical :: observed2(g_count, 1)
   real(dp) :: x2(g_count, 1)
   real(dp) :: prior1(2, 2)
   real(dp) :: prior_joint(3, 3)
   type(rng_state) :: rng
   type(jomo2hr_result) :: fit
   integer :: i
   integer :: g

   do g = 1, g_count
      x2(g, 1) = 1.0_dp
      y2(g, 1) = -0.3_dp + 0.12_dp * real(g - 3, dp) + 0.02_dp * cos(real(g, dp))
   end do
   do i = 1, n
      g = (i - 1) / per_group + 1
      cluster(i) = g
      x1(i, 1) = 1.0_dp
      x1(i, 2) = real(mod(i - 1, per_group), dp) / real(per_group - 1, dp)
      z(i, 1) = 1.0_dp
      y1(i, 1) = 0.5_dp + 0.4_dp * x1(i, 2) + 0.10_dp * y2(g, 1) + &
         (0.02_dp + 0.008_dp * real(g, dp)) * sin(real(i, dp))
      y1(i, 2) = -0.2_dp + 0.3_dp * x1(i, 2) - 0.08_dp * y2(g, 1) + &
         (0.03_dp + 0.006_dp * real(g, dp)) * cos(real(i, dp))
   end do
   observed1 = .true.
   observed2 = .true.
   observed1(5, 1) = .false.
   observed1(22, 2) = .false.
   observed1(33, :) = .false.
   observed2(4, 1) = .false.
   prior1 = 0.0_dp
   prior1(1, 1) = 1.0_dp
   prior1(2, 2) = 1.0_dp
   prior_joint = 0.0_dp
   prior_joint(1, 1) = 1.0_dp
   prior_joint(2, 2) = 1.0_dp
   prior_joint(3, 3) = 1.0_dp

   call rng_seed(rng, 16180339_i8)
   call jomo2conhr_mcmc(rng, y1, observed1, x1, z, cluster, y2, observed2, x2, 14, prior1, prior_joint, &
      5.0_dp, fit, hierarchy_df_prior=2.0_dp)

   if (.not. all(ieee_is_finite(fit%base%level1%continuous))) error stop "non-finite two-level heterogeneous imputation"
   if (.not. all(ieee_is_finite(fit%base%level2%continuous))) error stop "non-finite two-level heterogeneous level-2 state"
   if (.not. is_spd(fit%base%joint_covariance)) error stop "two-level heterogeneous joint covariance is not SPD"
   if (.not. is_spd(fit%hierarchy_scale)) error stop "two-level heterogeneous hierarchy scale is not SPD"
   if (.not. ieee_is_finite(fit%hierarchy_df) .or. fit%hierarchy_df < 2.0_dp) &
      error stop "two-level heterogeneous hierarchy degrees are invalid"
   do g = 1, g_count
      if (.not. is_spd(fit%cluster_covariance(g, :, :))) error stop "two-level cluster covariance is not SPD"
   end do
   if (maxval(abs(pack(fit%base%level1%continuous - y1, observed1))) > 1.0e-12_dp) &
      error stop "observed two-level heterogeneous level-1 values changed"
   if (maxval(abs(pack(fit%base%level2%continuous - y2, observed2))) > 1.0e-12_dp) &
      error stop "observed two-level heterogeneous level-2 values changed"

   print '(a)', "test_twolevel_heteroscedastic_continuous: PASS"
end program test_twolevel_heteroscedastic_continuous
