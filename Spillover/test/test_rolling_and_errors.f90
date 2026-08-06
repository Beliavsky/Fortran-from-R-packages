! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
program test_rolling_and_errors
   use iso_fortran_env, only : int64
   use spillover
   implicit none

   integer, parameter :: n = 240, k = 3
   real(dp) :: y(n, k), innovation(k), u
   real(dp), allocatable :: totals(:), nets(:, :)
   real(dp) :: ar(k, k, 1), sigma(k, k)
   type(var_model) :: model
   type(orthogonal_average_result) :: average
   integer(int64) :: state
   integer :: t, j, info
   character(len=200) :: message

   y = 0.0_dp
   state = 773_int64
   do t = 2, n
      do j = 1, k
         call next_uniform(state, u)
         innovation(j) = 0.30_dp * (u - 0.5_dp)
      end do
      y(t, 1) = 0.30_dp * y(t - 1, 1) + 0.20_dp * y(t - 1, 2) + innovation(1)
      y(t, 2) = 0.25_dp * y(t - 1, 2) + innovation(2)
      y(t, 3) = 0.20_dp * y(t - 1, 3) + 0.15_dp * y(t - 1, 1) + innovation(3)
   end do

   call rolling_total_spillover(y, 100, 1, 8, index_generalized, totals, info=info, message=message)
   call assert_true(info == spillover_success, trim(message))
   call rolling_net_spillover(y, 100, 1, 8, index_generalized, nets, info=info, message=message)
   call assert_true(info == spillover_success, trim(message))
   call assert_true(size(totals) == n - 99, 'rolling total length')
   call assert_true(size(nets, 1) == n - 99 .and. size(nets, 2) == k, 'rolling net shape')
   call assert_close(maxval(abs(sum(nets, dim=2))), 0.0_dp, 1.0e-9_dp, 'rolling net balance')

   call rolling_total_spillover(y, n + 1, 1, 8, index_generalized, totals, info=info, message=message)
   call assert_true(info == spillover_invalid_argument, 'invalid width status')

   ar = 0.0_dp
   sigma = 0.0_dp
   do j = 1, k
      sigma(j, j) = 1.0_dp
   end do
   call initialize_var_model(ar, sigma, model, info=info, message=message)
   call assert_true(info == spillover_success, trim(message))
   call orthogonal_average_exact(model, 2, average, exact_limit=2, info=info, message=message)
   call assert_true(info == spillover_iteration_limit, 'exact permutation dimension guard')

   print '(a)', 'test_rolling_and_errors: PASS'

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

end program test_rolling_and_errors
