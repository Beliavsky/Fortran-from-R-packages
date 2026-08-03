! SPDX-License-Identifier: MIT
program test_likelihood
   use invgamstochvol
   implicit none

   type(invgam_likelihood_result) :: fit
   real(dp) :: residuals(8)
   real(dp), parameter :: expected_total = -1.304456811364553_dp
   real(dp), parameter :: expected(0:7) = [ &
      -0.056352103210915_dp, 0.171303649959988_dp, &
      -0.524160952703820_dp, -0.157781518679341_dp, &
       0.224169473932321_dp, -0.716239494724050_dp, &
      -0.315184054102932_dp, 0.069788188164196_dp ]

   residuals = [0.20_dp, -0.10_dp, 0.35_dp, -0.25_dp, &
      0.05_dp, 0.40_dp, -0.30_dp, 0.15_dp]
   call lik_clo(residuals, 0.7_dp, 4.1_dp, 0.85_dp, fit, nit=20, niter=30)

   call check(fit%status == invgam_success, 'likelihood status')
   call check(abs(fit%total_loglik - expected_total) < 2.0e-12_dp, &
      'total log likelihood')
   call check(maxval(abs(fit%loglik - expected)) < 2.0e-12_dp, &
      'observation log likelihoods')
   call check(lbound(fit%all_st, 1) == 0, 'all_st lower bound')
   call check(ubound(fit%all_st, 1) == 8, 'all_st upper bound')
   call check(all(fit%all_st(1:8) > 0.0_dp), 'positive scale recursion')
   call check(size(fit%all_ctil, 1) == 8, 'coefficient row count')
   call check(size(fit%all_ctil, 2) == 21, 'coefficient column count')

   print '(a)', 'test_likelihood: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_likelihood
