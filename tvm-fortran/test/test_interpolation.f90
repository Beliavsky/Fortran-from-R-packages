! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
program test_interpolation
   use tvm_kinds, only : dp
   use tvm_interpolation, only : pchip_interpolator, build_pchip
   implicit none
   type(pchip_interpolator) :: interp
   real(dp) :: x(4), y(4), q(7), v(7)
   integer :: status, i

   x = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   y = [1.0_dp, 0.9_dp, 0.8_dp, 0.7_dp]
   call build_pchip(interp, x, y, status)
   if (status /= 0) error stop 1
   q = [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp]
   v = interp%evaluate_many(q)
   do i = 2, size(v)
      if (v(i) > v(i - 1)) error stop 1
   end do
   if (maxval(abs(interp%evaluate_many(x) - y)) > 1.0e-14_dp) error stop 1
   print '(a)', "test_interpolation: PASS"
end program test_interpolation
