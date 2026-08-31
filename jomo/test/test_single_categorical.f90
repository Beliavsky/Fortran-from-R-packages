program test_single_categorical
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_single_level, only : jomo1_result, jomo1cat_mcmc
   implicit none

   integer, parameter :: n = 24
   integer :: y(n, 2)
   logical :: observed(n, 2)
   integer :: n_levels(2)
   real(dp) :: x(n, 2)
   real(dp) :: prior(3, 3)
   type(rng_state) :: rng
   type(jomo1_result) :: fit
   integer :: i

   n_levels = [2, 3]
   do i = 1, n
      x(i, 1) = 1.0_dp
      x(i, 2) = (real(i, dp) - 12.5_dp) / 12.0_dp
      if (x(i, 2) > 0.0_dp) then
         y(i, 1) = 1
      else
         y(i, 1) = 2
      end if
      select case (mod(i, 3))
      case (0)
         y(i, 2) = 3
      case (1)
         y(i, 2) = 1
      case default
         y(i, 2) = 2
      end select
   end do
   observed = .true.
   observed(5, 1) = .false.
   observed(9, 2) = .false.
   observed(18, :) = .false.
   prior = 0.0_dp
   do i = 1, 3
      prior(i, i) = 1.0_dp
   end do

   call rng_seed(rng, 97531_i8)
   call jomo1cat_mcmc(rng, y, observed, n_levels, x, 20, prior, fit)

   if (.not. all(ieee_is_finite(fit%latent))) error stop "non-finite latent values"
   if (any(fit%categorical(:, 1) < 1) .or. any(fit%categorical(:, 1) > 2)) error stop "binary decode out of range"
   if (any(fit%categorical(:, 2) < 1) .or. any(fit%categorical(:, 2) > 3)) error stop "categorical decode out of range"
   do i = 1, n
      if (observed(i, 1)) then
         if (fit%categorical(i, 1) /= y(i, 1)) error stop "observed binary category changed"
      end if
      if (observed(i, 2)) then
         if (fit%categorical(i, 2) /= y(i, 2)) error stop "observed three-level category changed"
      end if
   end do
   if (abs(fit%omega(1, 1) - 1.0_dp) > 1.0e-12_dp) error stop "identified latent variance changed"
   if (abs(fit%omega(2, 2) - 1.0_dp) > 1.0e-12_dp) error stop "identified latent variance changed"
   if (abs(fit%omega(2, 3) - 0.5_dp) > 1.0e-12_dp) error stop "identified latent covariance changed"

   print '(a)', "test_single_categorical: PASS"
end program test_single_categorical
