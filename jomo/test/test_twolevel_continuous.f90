program test_twolevel_continuous
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_twolevel, only : jomo2_result, jomo2con_mcmc
   implicit none

   integer, parameter :: g_count = 8
   integer, parameter :: per_group = 5
   integer, parameter :: n = g_count * per_group
   real(dp) :: y1(n, 2)
   logical :: observed1(n, 2)
   real(dp) :: x1(n, 2)
   real(dp) :: z(n, 1)
   integer :: cluster(n)
   real(dp) :: y2(g_count, 1)
   logical :: observed2(g_count, 1)
   real(dp) :: x2(g_count, 2)
   real(dp) :: prior1(2, 2)
   real(dp) :: prior_joint(3, 3)
   type(rng_state) :: rng
   type(jomo2_result) :: fit
   integer :: i
   integer :: g
   real(dp) :: u1
   real(dp) :: u2

   do g = 1, g_count
      y2(g, 1) = -0.7_dp + 0.25_dp * real(g, dp) + 0.03_dp * sin(real(g, dp))
      x2(g, 1) = 1.0_dp
      x2(g, 2) = real(g - 1, dp) / real(g_count - 1, dp)
   end do
   do i = 1, n
      g = (i - 1) / per_group + 1
      cluster(i) = g
      x1(i, 1) = 1.0_dp
      x1(i, 2) = real(mod(i - 1, per_group), dp) / real(per_group - 1, dp)
      z(i, 1) = 1.0_dp
      u1 = 0.18_dp * y2(g, 1)
      u2 = -0.12_dp * y2(g, 1)
      y1(i, 1) = 1.0_dp + 0.6_dp * x1(i, 2) + u1 + 0.04_dp * sin(real(i, dp))
      y1(i, 2) = -0.3_dp + 0.2_dp * x1(i, 2) + u2 + 0.04_dp * cos(real(i, dp))
   end do
   observed1 = .true.
   observed2 = .true.
   observed1(6, 1) = .false.
   observed1(31, 2) = .false.
   observed2(4, 1) = .false.
   prior1 = 0.0_dp
   prior_joint = 0.0_dp
   prior1(1, 1) = 1.0_dp
   prior1(2, 2) = 1.0_dp
   do i = 1, 3
      prior_joint(i, i) = 1.0_dp
   end do

   call rng_seed(rng, 87654321_i8)
   call jomo2con_mcmc(rng, y1, observed1, x1, z, cluster, y2, observed2, x2, 18, prior1, prior_joint, fit)

   if (.not. all(ieee_is_finite(fit%level1%continuous))) error stop "non-finite level-1 imputation"
   if (.not. all(ieee_is_finite(fit%level2%continuous))) error stop "non-finite level-2 imputation"
   if (.not. is_spd(fit%level1%omega)) error stop "level-1 covariance is not SPD"
   if (.not. is_spd(fit%joint_covariance)) error stop "joint covariance is not SPD"
   if (maxval(abs(pack(fit%level1%continuous - y1, observed1))) > 1.0e-12_dp) error stop "observed level-1 changed"
   if (maxval(abs(pack(fit%level2%continuous - y2, observed2))) > 1.0e-12_dp) error stop "observed level-2 changed"

   print '(a)', "test_twolevel_continuous: PASS"
end program test_twolevel_continuous
