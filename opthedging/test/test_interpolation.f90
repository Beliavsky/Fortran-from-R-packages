! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
program test_interpolation
   use opthedging, only : dp, interpolation1d
   implicit none

   real(dp) :: identity(10)
   real(dp) :: constant(1)
   integer :: i

   identity = [(real(i, dp), i=1, 10)]
   call check_close(interpolation1d(2.45_dp, identity, 1.0_dp, 10.0_dp), &
      2.45_dp, 1.0e-14_dp, "interior interpolation")
   call check_close(interpolation1d(0.25_dp, identity, 1.0_dp, 10.0_dp), &
      0.25_dp, 1.0e-14_dp, "left extrapolation")
   call check_close(interpolation1d(11.75_dp, identity, 1.0_dp, 10.0_dp), &
      11.75_dp, 1.0e-14_dp, "right extrapolation")
   constant = 7.0_dp
   call check_close(interpolation1d(100.0_dp, constant, 0.0_dp, 1.0_dp), &
      7.0_dp, 0.0_dp, "single point")

   print '(a)', 'test_interpolation: PASS'

contains

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual
      real(dp), intent(in) :: expected
      real(dp), intent(in) :: tolerance
      character(len=*), intent(in) :: label

      if (abs(actual - expected) > tolerance) then
         print '(a)', 'failed: ' // label
         print '(a,es24.16)', 'actual:   ', actual
         print '(a,es24.16)', 'expected: ', expected
         error stop 1
      end if
   end subroutine check_close

end program test_interpolation
