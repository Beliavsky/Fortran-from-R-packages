! SPDX-License-Identifier: GPL-2.0-or-later
program test_dvi_rounding
   use ttr, only: dp, dvi
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   real(dp) :: x(512), expected(512, 3), actual(512, 3)
   integer :: unit, i, j
   open(newunit = unit, file = 'test/dvi_reference.dat', status = 'old', action = 'read')
   do i = 1, size(x)
      read(unit, *) x(i), expected(i, :)
   end do
   close(unit)
   actual = dvi(x, 21, weights = [0.3_dp, 0.7_dp], smooth = 4, &
      magnitude = [7, 19, 4], stretch = [5, 13, 3], exact_multiplier = 0.3_dp)
   do j = 1, 3
      do i = 1, size(x)
         if (ieee_is_nan(expected(i, j))) then
            if (ieee_is_nan(actual(i, j))) cycle
         else if (.not. ieee_is_nan(actual(i, j))) then
            if (abs(actual(i, j) - expected(i, j)) < 1.0e-11_dp) cycle
         end if
         print *, 'DVI recurrence regression row/column ', i, j, actual(i, j), expected(i, j)
         error stop 'DVI percent-rank tie regression'
      end do
   end do
   print *, 'DVI sliding-sum/rank regression: 1536 genuine R values passed'
end program test_dvi_rounding
