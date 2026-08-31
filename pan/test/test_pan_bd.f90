! SPDX-License-Identifier: GPL-3.0-only
program test_pan_bd
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use pan, only : PAN_OK, dp, pan_bd_mcmc, pan_bd_prior, pan_bd_result
   use pan_linalg, only : is_spd
   implicit none

   type(pan_bd_prior) :: prior
   type(pan_bd_result) :: fit
   integer :: a
   integer :: i
   integer :: subj(16)
   integer :: xcol(2)
   integer :: zcol(1)
   real(dp) :: pred(16, 2)
   real(dp) :: y(16, 2)
   real(dp) :: nan

   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   subj = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
   xcol = [1, 2]
   zcol = [1]

   do i = 1, 16
      pred(i, 1) = 1.0_dp
      pred(i, 2) = real(mod(i - 1, 4), dp)
      y(i, 1) = 2.0_dp + 0.2_dp * pred(i, 2) + 0.03_dp * real(i, dp)
      y(i, 2) = -1.0_dp + 0.5_dp * pred(i, 2) - 0.02_dp * real(i, dp)
   end do
   y(4, 1) = nan
   y(10, 2) = nan

   prior%a = 2.0_dp
   allocate(prior%binv(2, 2), prior%c(2), prior%dinv(1, 1, 2))
   prior%binv = reshape([1.0_dp, 0.1_dp, 0.1_dp, 1.0_dp], [2, 2])
   prior%c = [1.0_dp, 1.0_dp]
   prior%dinv(1, 1, :) = [1.0_dp, 1.5_dp]

   call pan_bd_mcmc(y, subj, pred, xcol, zcol, prior, 8888, 5, fit)

   call assert_true(fit%status == PAN_OK, "pan.bd sampler status")
   call assert_true(all(shape(fit%psi) == [1, 1, 2, 5]), "pan.bd psi chain shape")
   call assert_true(.not. any(ieee_is_nan(fit%y)), "pan.bd imputes all missing values")
   call assert_true(is_spd(fit%last%sigma), "pan.bd residual covariance SPD")
   do a = 1, 2
      call assert_true(is_spd(fit%last%psi(:, :, a)), "pan.bd random-effect block SPD")
   end do

   call finish_tests()

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Test predicate expected to be true.
      character(len=*), intent(in) :: message !! Diagnostic printed if the predicate is false.

      if (.not. condition) then
         write (*, '(a)') "FAIL: " // trim(message)
         error stop 1
      end if
   end subroutine assert_true


   subroutine finish_tests()
   end subroutine finish_tests

end program test_pan_bd
