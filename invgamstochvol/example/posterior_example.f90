! SPDX-License-Identifier: MIT
program posterior_example
   use, intrinsic :: iso_fortran_env, only : int64
   use invgamstochvol
   implicit none

   type(invgam_likelihood_result) :: fit
   real(dp) :: residuals(8)
   real(dp), allocatable :: inverse_volatility(:)
   integer :: t, status

   residuals = [0.20_dp, -0.10_dp, 0.35_dp, -0.25_dp, &
      0.05_dp, 0.40_dp, -0.30_dp, 0.15_dp]
   call lik_clo(residuals, 0.7_dp, 4.1_dp, 0.85_dp, fit, nit=30, niter=40)
   call draw_k0(fit, 4.1_dp, 0.85_dp, 0.7_dp, inverse_volatility, &
      seed=12345_int64, status=status)

   write (*, '(a,i0)') 'draw status = ', status
   write (*, '(a)') 't       inverse volatility       volatility / b2'
   do t = 1, size(inverse_volatility)
      write (*, '(i0,3x,f18.10,3x,f18.10)') t, inverse_volatility(t), &
         1.0_dp / (0.7_dp * inverse_volatility(t))
   end do
end program posterior_example
