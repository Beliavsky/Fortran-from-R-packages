! SPDX-License-Identifier: GPL-3.0-only
program pan_example
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use pan, only : PAN_OK, dp, ecme_fit, ecme_result, pan_mcmc, pan_prior, pan_result
   implicit none

   type(ecme_result) :: ml
   type(pan_prior) :: prior
   type(pan_result) :: imp
   integer :: i
   integer :: occ(12)
   integer :: subj(12)
   integer :: xcol(2)
   integer :: zcol(1)
   real(dp) :: pred(12, 2)
   real(dp) :: y_complete(12)
   real(dp) :: y_missing(12, 1)
   real(dp) :: nan

   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   subj = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
   occ = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]
   xcol = [1, 2]
   zcol = [1]

   do i = 1, 12
      pred(i, 1) = 1.0_dp
      pred(i, 2) = real(occ(i) - 1, dp)
      y_complete(i) = 1.0_dp + 0.4_dp * pred(i, 2) + 0.08_dp * real(subj(i) - 2, dp)
   end do

   y_missing(:, 1) = y_complete
   y_missing(5, 1) = nan
   y_missing(11, 1) = nan

   prior%a = 1.0_dp
   prior%binv = reshape([1.0_dp], [1, 1])
   prior%c = 1.0_dp
   prior%dinv = reshape([1.0_dp], [1, 1])

   call pan_mcmc(y_missing, subj, pred, xcol, zcol, prior, 13579, 20, imp)
   if (imp%status /= PAN_OK) error stop imp%message

   write (*, '(a,2f10.4)') "imputed rows 5 and 11:", imp%y(5, 1), imp%y(11, 1)
   write (*, '(a,f10.4)') "final residual variance:", imp%last%sigma(1, 1)
   write (*, '(a,f10.4)') "final random-intercept variance:", imp%last%psi(1, 1)

   call ecme_fit(y_complete, subj, occ, pred, xcol, ml, zcol=zcol)
   if (ml%status /= PAN_OK) error stop ml%message

   write (*, '(a,2f10.4)') "ecme fixed effects:", ml%beta
   write (*, '(a,f10.4)') "ecme residual variance:", ml%sigma2
   write (*, '(a,f10.4)') "ecme random-intercept variance:", ml%psi(1, 1)
end program pan_example
