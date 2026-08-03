! SPDX-License-Identifier: MIT
program test_draw
   use, intrinsic :: iso_fortran_env, only : int64
   use invgamstochvol
   implicit none

   type(invgam_likelihood_result) :: fit
   real(dp) :: residuals(10)
   real(dp), allocatable :: draw1(:), draw2(:), draw3(:)
   integer :: status

   residuals = [0.15_dp, -0.20_dp, 0.05_dp, 0.30_dp, -0.10_dp, &
      0.12_dp, -0.18_dp, 0.22_dp, 0.08_dp, -0.05_dp]
   call lik_clo(residuals, 0.8_dp, 5.0_dp, 0.7_dp, fit, nit=25, niter=35)
   call check(fit%status == invgam_success, 'fit before draw')

   call draw_k0(fit, 5.0_dp, 0.7_dp, 0.8_dp, draw1, &
      seed=24681357_int64, status=status)
   call check(status == invgam_success, 'first posterior draw')
   call draw_k0(fit, 5.0_dp, 0.7_dp, 0.8_dp, draw2, &
      seed=24681357_int64, status=status)
   call check(status == invgam_success, 'repeated posterior draw')
   call draw_k0(fit, 5.0_dp, 0.7_dp, 0.8_dp, draw3, &
      seed=97531_int64, status=status)

   call check(size(draw1) == size(residuals), 'posterior draw length')
   call check(all(draw1 > 0.0_dp), 'positive inverse volatilities')
   call check(maxval(abs(draw1 - draw2)) < tiny(1.0_dp), &
      'fixed seed reproducibility')
   call check(maxval(abs(draw1 - draw3)) > 1.0e-8_dp, &
      'different seeds change draw')

   print '(a)', 'test_draw: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_draw
