! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the computational code of R package fda 6.3.0.
program smooth_sine
   use fda, only : dp, basis_type, make_bspline_basis, smooth_basis, smooth_result_type
   implicit none
   type(basis_type) :: basis
   type(smooth_result_type) :: fit
   real(dp) :: breaks(6), t(21), y(21, 1)
   integer :: i, info

   breaks = [0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp, 1.0_dp]
   do i = 1, size(t)
      t(i) = real(i - 1, dp) / real(size(t) - 1, dp)
      y(i, 1) = sin(2.0_dp * acos(-1.0_dp) * t(i))
   end do
   call make_bspline_basis(breaks, 4, basis, info)
   if (info /= 0) error stop "failed to create B-spline basis"
   call smooth_basis(t, y, basis, 1.0e-5_dp, 2, fit, info)
   if (info /= 0) error stop "smoothing failed"
   write (*, '(a,f10.5)') "effective df: ", fit%df
   write (*, '(a,es12.4)') "SSE:          ", fit%sse
end program smooth_sine
