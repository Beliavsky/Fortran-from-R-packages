program test_twolevel_categorical
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_twolevel, only : jomo2_result, jomo2_mixed_mcmc
   implicit none

   integer, parameter :: g_count = 10
   integer, parameter :: per_group = 4
   integer, parameter :: n = g_count * per_group
   real(dp) :: y1(n, 1)
   logical :: observed1(n, 1)
   integer :: empty_cat1(n, 0)
   logical :: empty_cat1_obs(n, 0)
   integer :: levels1(0)
   real(dp) :: x1(n, 2)
   real(dp) :: z(n, 1)
   integer :: cluster(n)
   real(dp) :: empty_con2(g_count, 0)
   logical :: empty_con2_obs(g_count, 0)
   integer :: y2_cat(g_count, 1)
   logical :: observed2_cat(g_count, 1)
   integer :: levels2(1)
   real(dp) :: x2(g_count, 1)
   real(dp) :: prior1(1, 1)
   real(dp) :: prior_joint(2, 2)
   type(rng_state) :: rng
   type(jomo2_result) :: fit
   integer :: i
   integer :: g

   levels2 = 2
   do g = 1, g_count
      y2_cat(g, 1) = merge(1, 2, mod(g, 2) == 0)
      x2(g, 1) = 1.0_dp
   end do
   do i = 1, n
      g = (i - 1) / per_group + 1
      cluster(i) = g
      x1(i, 1) = 1.0_dp
      x1(i, 2) = real(mod(i - 1, per_group), dp) / real(per_group - 1, dp)
      z(i, 1) = 1.0_dp
      y1(i, 1) = 0.4_dp + 0.3_dp * x1(i, 2) + merge(0.12_dp, -0.12_dp, y2_cat(g, 1) == 1) + &
         0.03_dp * sin(real(i, dp))
   end do
   observed1 = .true.
   observed1(9, 1) = .false.
   observed2_cat = .true.
   observed2_cat(7, 1) = .false.
   prior1(1, 1) = 1.0_dp
   prior_joint = 0.0_dp
   prior_joint(1, 1) = 1.0_dp
   prior_joint(2, 2) = 1.0_dp

   call rng_seed(rng, 99887766_i8)
   call jomo2_mixed_mcmc(rng, y1, observed1, empty_cat1, empty_cat1_obs, levels1, x1, z, cluster, &
      empty_con2, empty_con2_obs, y2_cat, observed2_cat, levels2, x2, 14, prior1, prior_joint, fit)

   if (.not. all(ieee_is_finite(fit%level1%continuous))) error stop "non-finite mixed two-level imputation"
   if (.not. is_spd(fit%joint_covariance)) error stop "mixed two-level joint covariance is not SPD"
   if (abs(fit%joint_covariance(2, 2) - 1.0_dp) > 1.0e-12_dp) error stop "binary latent variance identification changed"
   do g = 1, g_count
      if (observed2_cat(g, 1)) then
         if (fit%level2%categorical(g, 1) /= y2_cat(g, 1)) error stop "observed level-2 category changed"
      end if
   end do
   if (fit%level2%categorical(7, 1) < 1 .or. fit%level2%categorical(7, 1) > 2) error stop "invalid imputed category"

   print '(a)', "test_twolevel_categorical: PASS"
end program test_twolevel_categorical
