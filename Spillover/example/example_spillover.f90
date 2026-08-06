! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
program example_spillover
   use iso_fortran_env, only : int64
   use spillover
   implicit none

   integer, parameter :: n = 600, k = 3
   real(dp) :: y(n, k), innovation(k), u
   type(var_model) :: model
   type(spillover_result) :: result
   integer(int64) :: state
   integer :: t, j, info
   character(len=200) :: message

   y = 0.0_dp
   state = 20260803_int64
   do t = 2, n
      do j = 1, k
         state = modulo(16807_int64 * state, 2147483647_int64)
         u = real(state, dp) / 2147483647.0_dp
         innovation(j) = 0.35_dp * (u - 0.5_dp)
      end do
      y(t, 1) = 0.35_dp * y(t - 1, 1) + 0.30_dp * y(t - 1, 2) + innovation(1)
      y(t, 2) = 0.25_dp * y(t - 1, 2) + 0.15_dp * y(t - 1, 3) + innovation(2)
      y(t, 3) = 0.20_dp * y(t - 1, 3) + 0.10_dp * y(t - 1, 1) + innovation(3)
   end do

   call fit_var(y, 2, model, var_const, info, message)
   if (info /= spillover_success) error stop trim(message)
   call g_spillover(model, 10, .true., result, info, message)
   if (info /= spillover_success) error stop trim(message)

   write(*, '(a,f10.4)') 'Generalized total connectedness: ', result%total
   write(*, '(a)') 'Variable       FROM          TO         NET'
   do j = 1, k
      write(*, '(i5,3f12.4)') j, result%from(j), result%to(j), result%net(j)
   end do
end program example_spillover
