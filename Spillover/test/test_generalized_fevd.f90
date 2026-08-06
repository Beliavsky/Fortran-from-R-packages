! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
program test_generalized_fevd
   use spillover
   implicit none

   real(dp) :: ar(2, 2, 1), sigma(2, 2)
   real(dp), allocatable :: fevd(:, :, :), raw(:, :, :)
   type(var_model) :: model
   type(spillover_result) :: result
   integer :: info, h
   character(len=200) :: message

   ar = 0.0_dp
   sigma = reshape([1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
   call initialize_var_model(ar, sigma, model, info=info, message=message)
   call assert_true(info == spillover_success, trim(message))

   call generalized_fevd(model, 3, .true., fevd, info, message)
   call assert_true(info == spillover_success, trim(message))
   do h = 1, 3
      call assert_close(fevd(1, 1, h), 0.8_dp, 1.0e-13_dp, 'normalized own share')
      call assert_close(fevd(1, 2, h), 0.2_dp, 1.0e-13_dp, 'normalized cross share')
      call assert_close(sum(fevd(1, :, h)), 1.0_dp, 1.0e-13_dp, 'row normalization')
   end do

   call generalized_fevd(model, 1, .false., raw, info, message)
   call assert_true(info == spillover_success, trim(message))
   call assert_close(raw(1, 1, 1), 1.0_dp, 1.0e-13_dp, 'raw own share')
   call assert_close(raw(1, 2, 1), 0.25_dp, 1.0e-13_dp, 'raw cross share')

   call generalized_spillover(model, 1, .true., result, info, message)
   call assert_true(info == spillover_success, trim(message))
   call assert_close(result%total, 20.0_dp, 1.0e-12_dp, 'generalized total index')
   call assert_close(maxval(abs(result%net)), 0.0_dp, 1.0e-12_dp, 'symmetric net index')

   print '(a)', 'test_generalized_fevd: PASS'

contains

   subroutine assert_true(condition, text)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: text
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // text
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, text)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: text
      if (abs(actual - expected) > tolerance) then
         write(*, '(a,2(1x,es14.6))') 'FAIL: ' // text, actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_generalized_fevd
