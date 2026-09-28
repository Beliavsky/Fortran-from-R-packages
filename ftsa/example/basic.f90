! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
program basic
   use ftsa, only : dp, ftsm_fit, ftsa_forecast_ftsm, ftsa_ftsm
   implicit none
   type(ftsm_fit) :: fit
   real(dp) :: y(3, 5)
   real(dp), allocatable :: forecast(:, :)
   integer :: j, status

   do j = 1, 5
      y(:, j) = [1.0_dp, 2.0_dp, 3.0_dp] + real(j - 3, dp) * [0.2_dp, -0.1_dp, 0.3_dp]
   end do
   call ftsa_ftsm(y, 1, fit, status)
   if (status /= 0) error stop 'ftsa_ftsm failed'
   call ftsa_forecast_ftsm(fit, 2, forecast, status)
   if (status /= 0) error stop 'ftsa_forecast_ftsm failed'
   print '(a,3f10.4)', 'One-step functional forecast:', forecast(:, 1)
end program basic
