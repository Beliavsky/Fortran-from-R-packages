program test_single_continuous
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_single_level, only : jomo1_result, jomo1con_mcmc
   implicit none

   integer, parameter :: n = 16
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
      y(i, 1) = 1.0_dp + 0.8_dp * x(i, 2) + 0.15_dp * sin(real(i, dp))
      y(i, 2) = -0.5_dp + 0.4_dp * x(i, 2) + 0.10_dp * cos(real(i, dp))
   end do
   observed = .true.
   observed(3, 1) = .false.
   observed(7, 2) = .false.
   observed(12, :) = .false.
   prior = 0.0_dp
   prior(1, 1) = 1.0_dp
   prior(2, 2) = 1.0_dp

   call rng_seed(rng, 24681357_i8)
   call jomo1con_mcmc(rng, y, observed, x, 30, prior, fit, store_chain=.true.)

   if (fit%iterations /= 30) error stop "wrong iteration count"
   if (.not. all(ieee_is_finite(fit%continuous))) error stop "non-finite continuous imputation"
   if (.not. all(ieee_is_finite(fit%beta))) error stop "non-finite beta"
   if (.not. is_spd(fit%omega)) error stop "omega is not SPD"
   if (maxval(abs(pack(fit%continuous - y, observed))) > 1.0e-12_dp) error stop "observed values changed"
   if (size(fit%beta_chain, 3) /= 30) error stop "beta chain was not stored"
   if (abs(fit%continuous(3, 1)) > 100.0_dp) error stop "implausible missing-value draw"

   print '(a)', "test_single_continuous: PASS"
end program test_single_continuous
