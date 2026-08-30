! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
program fit_var_example
   use r_kinds, only : dp
   use vars, only : var_model, forecast_result, fit_var, forecast_var, var_const, vars_success
   implicit none

   type(var_model) :: model
   type(forecast_result) :: forecast
   real(dp) :: y(12, 2)
   integer :: i, info

   y(1, :) = [0.2_dp, -0.1_dp]
   do i = 2, size(y, 1)
      y(i, 1) = 0.7_dp + 0.45_dp * y(i - 1, 1) + 0.10_dp * y(i - 1, 2) &
         + 0.02_dp * real(mod(i, 3) - 1, dp)
      y(i, 2) = -0.2_dp - 0.15_dp * y(i - 1, 1) + 0.35_dp * y(i - 1, 2) &
         + 0.01_dp * real(mod(i + 1, 4) - 2, dp)
   end do

   call fit_var(y, 1, var_const, model, info)
   if (info /= vars_success) error stop "fit_var failed"

   call forecast_var(model, 3, 0.95_dp, forecast, info)
   if (info /= vars_success) error stop "forecast_var failed"

   write (*, '(a)') "VAR(1) coefficient matrix:"
   do i = 1, model%k
      write (*, '(3f14.7)') model%coef(i, :)
   end do
   write (*, '(a)') "Three-step point forecasts:"
   do i = 1, 3
      write (*, '(2f14.7)') forecast%point(i, :)
   end do
end program fit_var_example
