! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
program test_var
   use iso_fortran_env, only : int64
   use spillover
   implicit none

   integer, parameter :: n = 800, k = 2
   real(dp) :: y(n, k), a(k, k), c(k), e(k), u
   real(dp), allocatable :: phi(:, :, :)
   type(var_model) :: model
   integer(int64) :: state
   integer :: t, info
   character(len=200) :: message

   a = reshape([0.45_dp, -0.10_dp, 0.20_dp, 0.35_dp], [k, k])
   c = [0.10_dp, -0.05_dp]
   y = 0.0_dp
   state = 12457_int64
   do t = 2, n
      call next_uniform(state, u)
      e(1) = 0.30_dp * (u - 0.5_dp)
      call next_uniform(state, u)
      e(2) = 0.25_dp * (u - 0.5_dp)
      y(t, :) = c + matmul(a, y(t - 1, :)) + e
   end do

   call fit_var(y, 1, model, var_const, info, message)
   call assert_true(info == spillover_success, trim(message))
   call assert_close(maxval(abs(model%ar(:, :, 1) - a)), 0.0_dp, 0.060_dp, &
      'VAR coefficient recovery')
   call assert_close(maxval(abs(model%intercept - c)), 0.0_dp, 0.015_dp, &
      'VAR intercept recovery')
   call assert_true(all([(model%sigma(t, t) > 0.0_dp, t = 1, k)]), &
      'positive residual variances')

   call ma_coefficients(model, 4, phi, info, message)
   call assert_true(info == spillover_success, trim(message))
   call assert_close(maxval(abs(phi(:, :, 1) - reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [k, k]))), &
      0.0_dp, 1.0e-13_dp, 'Phi zero')
   call assert_close(maxval(abs(phi(:, :, 2) - model%ar(:, :, 1))), 0.0_dp, 1.0e-13_dp, 'Phi one')
   call assert_close(maxval(abs(phi(:, :, 3) - matmul(model%ar(:, :, 1), model%ar(:, :, 1)))), &
      0.0_dp, 1.0e-12_dp, 'Phi two')

   print '(a)', 'test_var: PASS'

contains

   subroutine next_uniform(state, value)
      integer(int64), intent(inout) :: state
      real(dp), intent(out) :: value
      state = modulo(16807_int64 * state, 2147483647_int64)
      value = real(state, dp) / 2147483647.0_dp
   end subroutine next_uniform

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

end program test_var
