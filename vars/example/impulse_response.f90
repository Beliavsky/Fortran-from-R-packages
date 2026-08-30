! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
program impulse_response_example
   use r_kinds, only : dp
   use vars, only : var_model, fit_var, impulse_response, fevd_var, var_const, vars_success
   implicit none

   type(var_model) :: model
   real(dp) :: y(16, 2)
   real(dp), allocatable :: irf(:, :, :), fevd(:, :, :)
   integer :: i, h, info

   y(1, :) = [0.1_dp, -0.2_dp]
   do i = 2, size(y, 1)
      y(i, 1) = 0.5_dp + 0.50_dp * y(i - 1, 1) + 0.08_dp * y(i - 1, 2) &
         + 0.015_dp * real(mod(i, 4) - 2, dp)
      y(i, 2) = -0.3_dp - 0.12_dp * y(i - 1, 1) + 0.42_dp * y(i - 1, 2) &
         + 0.012_dp * real(mod(i, 5) - 2, dp)
   end do

   call fit_var(y, 1, var_const, model, info)
   if (info /= vars_success) error stop "fit_var failed"

   call impulse_response(model, 5, .true., .false., irf, info)
   if (info /= vars_success) error stop "impulse_response failed"
   call fevd_var(model, 5, fevd, info)
   if (info /= vars_success) error stop "fevd_var failed"

   write (*, '(a)') "Orthogonalized response of variable 1 to shock 1:"
   do h = 1, size(irf, 3)
      write (*, '(i3,2x,f14.7)') h - 1, irf(1, 1, h)
   end do
   write (*, '(a,f12.7)') "Five-step FEVD share from shock 2: ", fevd(1, 2, 5)
end program impulse_response_example
