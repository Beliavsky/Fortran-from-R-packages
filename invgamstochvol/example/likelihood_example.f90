! SPDX-License-Identifier: MIT
program likelihood_example
   use invgamstochvol
   implicit none

   type(invgam_likelihood_result) :: fit
   real(dp) :: residuals(8)
   integer :: t

   residuals = [0.20_dp, -0.10_dp, 0.35_dp, -0.25_dp, &
      0.05_dp, 0.40_dp, -0.30_dp, 0.15_dp]
   call lik_clo(residuals, 0.7_dp, 4.1_dp, 0.85_dp, fit, nit=20, niter=30)

   write (*, '(a,f18.10)') 'total log likelihood = ', fit%total_loglik
   write (*, '(a)') 't       contribution'
   do t = 0, size(fit%loglik) - 1
      write (*, '(i0,3x,f18.10)') t + 1, fit%loglik(t)
   end do
end program likelihood_example
