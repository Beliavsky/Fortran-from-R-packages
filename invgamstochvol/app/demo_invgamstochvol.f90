! SPDX-License-Identifier: MIT
program demo_invgamstochvol
   use, intrinsic :: iso_fortran_env, only : int64
   use invgamstochvol
   implicit none

   type(invgam_likelihood_result) :: fit
   real(dp), parameter :: b2 = 0.7_dp, nu = 4.1_dp, rho = 0.85_dp
   real(dp) :: residuals(8)
   real(dp), allocatable :: inverse_volatility(:)
   integer :: status

   residuals = [0.20_dp, -0.10_dp, 0.35_dp, -0.25_dp, &
      0.05_dp, 0.40_dp, -0.30_dp, 0.15_dp]

   call lik_clo(residuals, b2, nu, rho, fit, nit=30, niter=40)
   if (fit%status /= invgam_success) then
      write (*, '(a)') trim(fit%message)
      error stop 1
   end if

   call draw_k0(fit, nu, rho, b2, inverse_volatility, &
      seed=20260731_int64, status=status)
   if (status /= invgam_success) error stop 2

   write (*, '(a,f16.8)') 'Exact log likelihood: ', fit%total_loglik
   write (*, '(a,f16.8)') 'Mean posterior volatility draw: ', &
      sum(1.0_dp / (b2 * inverse_volatility)) / real(size(inverse_volatility), dp)
end program demo_invgamstochvol
