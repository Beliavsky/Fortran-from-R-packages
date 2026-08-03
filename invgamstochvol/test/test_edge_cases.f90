! SPDX-License-Identifier: MIT
program test_edge_cases
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use invgamstochvol
   implicit none

   type(invgam_likelihood_result) :: fit
   real(dp) :: short_series(2), residuals(5), bad(5)
   real(dp), allocatable :: draw(:)
   integer :: status

   short_series = [0.1_dp, -0.2_dp]
   call lik_clo(short_series, 1.0_dp, 4.0_dp, 0.8_dp, fit)
   call check(fit%status == invgam_invalid_argument, 'short series rejected')

   residuals = [0.1_dp, -0.1_dp, 0.2_dp, -0.2_dp, 0.05_dp]
   call lik_clo(residuals, 0.9_dp, 4.5_dp, 0.0_dp, fit, nit=0, niter=10)
   call check(fit%status == invgam_success, 'zero persistence and nit zero')
   call check(abs(fit%total_loglik) < huge(1.0_dp), 'finite edge likelihood')
   call draw_k0(fit, 4.5_dp, 0.0_dp, 0.9_dp, draw, status=status)
   call check(status == invgam_success, 'zero-persistence posterior draw')
   call check(all(draw > 0.0_dp), 'zero-persistence draw positive')

   bad = residuals
   bad(3) = ieee_value(0.0_dp, ieee_quiet_nan)
   call lik_clo(bad, 0.9_dp, 4.5_dp, 0.5_dp, fit)
   call check(fit%status == invgam_nonfinite_input, 'NaN rejected')

   call lik_clo(residuals, -1.0_dp, 4.5_dp, 0.5_dp, fit)
   call check(fit%status == invgam_invalid_argument, 'negative b2 rejected')

   print '(a)', 'test_edge_cases: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_edge_cases
