program test_multilevel_categorical
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_multilevel, only : jomo1ran_result, jomo1ran_mixed_mcmc
   implicit none

   integer, parameter :: g_count = 6
   integer, parameter :: per_group = 6
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
   type(jomo1ran_result) :: fit
   integer :: i
   integer :: g

   n_levels = 2
   do i = 1, n
      g = (i - 1) / per_group + 1
      cluster(i) = g
      x(i, 1) = 1.0_dp
      x(i, 2) = -1.0_dp + 2.0_dp * real(mod(i - 1, per_group), dp) / real(per_group - 1, dp)
      z(i, 1) = 1.0_dp
      y_con(i, 1) = 0.4_dp + 0.3_dp * x(i, 2) + 0.08_dp * real(g - 3, dp) + 0.03_dp * sin(real(i, dp))
      if (x(i, 2) + 0.20_dp * real(g - 3, dp) > 0.0_dp) then
         y_cat(i, 1) = 1
      else
         y_cat(i, 1) = 2
      end if
   end do
   con_observed = .true.
   cat_observed = .true.
   con_observed(8, 1) = .false.
   cat_observed(17, 1) = .false.
   con_observed(31, 1) = .false.
   cat_observed(31, 1) = .false.
   prior = 0.0_dp
   prior_u = 0.0_dp
   prior(1, 1) = 1.0_dp
   prior(2, 2) = 1.0_dp
   prior_u(1, 1) = 1.0_dp
   prior_u(2, 2) = 1.0_dp

   call rng_seed(rng, 31415926_i8)
   call jomo1ran_mixed_mcmc(rng, y_con, con_observed, y_cat, cat_observed, n_levels, x, z, cluster, 18, &
      prior, prior_u, fit)

   if (.not. all(ieee_is_finite(fit%level1%latent))) error stop "multilevel categorical latent values are not finite"
   if (.not. is_spd(fit%random_covariance)) error stop "multilevel categorical random covariance is not SPD"
   if (abs(fit%level1%omega(2, 2) - 1.0_dp) > 1.0e-12_dp) error stop "binary latent variance changed"
   do i = 1, n
      if (cat_observed(i, 1)) then
         if (fit%level1%categorical(i, 1) /= y_cat(i, 1)) error stop "observed multilevel category changed"
      end if
   end do
   if (fit%level1%categorical(17, 1) < 1 .or. fit%level1%categorical(17, 1) > 2) &
      error stop "multilevel categorical imputation out of range"

   print '(a)', "test_multilevel_categorical: PASS"
end program test_multilevel_categorical
