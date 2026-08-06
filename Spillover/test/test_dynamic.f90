! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
program test_dynamic
   use iso_fortran_env, only : int64
   use spillover
   implicit none

   integer, parameter :: n = 420, k = 2
   real(dp) :: y(n, k), a(k, k), innovation(k), u
   real(dp) :: early_mean, late_mean
   type(dynamic_spillover_result) :: result
   integer(int64) :: state
   integer :: t, info
   character(len=200) :: message

   y = 0.0_dp
   state = 9137_int64
   do t = 2, n
      if (t <= 210) then
         a = reshape([0.25_dp, 0.0_dp, 0.0_dp, 0.20_dp], [k, k])
      else
         a = reshape([0.25_dp, 0.0_dp, 0.65_dp, 0.20_dp], [k, k])
      end if
      call next_uniform(state, u)
      innovation(1) = 0.40_dp * (u - 0.5_dp)
      call next_uniform(state, u)
      innovation(2) = 0.40_dp * (u - 0.5_dp)
      y(t, :) = matmul(a, y(t - 1, :)) + innovation
   end do

   call dynamic_spillover(y, 100, 1, 10, result, info=info, message=message)
   call assert_true(info == spillover_success, trim(message))
   call assert_true(result%n_windows == n - 99, 'rolling window count')
   call assert_true(result%n_pairs == 1, 'pair count')
   call assert_close(maxval(abs(sum(result%net, dim=2))), 0.0_dp, 1.0e-9_dp, &
      'net spillovers sum to zero')

   early_mean = sum(result%total(1:80)) / 80.0_dp
   late_mean = sum(result%total(result%n_windows - 79:result%n_windows)) / 80.0_dp
   call assert_true(late_mean > early_mean + 5.0_dp, 'dynamic index detects cross-lag regime')

   print '(a,2(1x,f10.4))', 'test_dynamic: PASS; early/late totals', early_mean, late_mean

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

end program test_dynamic
