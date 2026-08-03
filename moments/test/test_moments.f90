! SPDX-License-Identifier: GPL-2.0-or-later
program test_moments
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use moments, only : dp, moment, all_moments, skewness, kurtosis, geary
   implicit none

   real(dp) :: x(10), xnan(5), m(5, 2)
   real(dp), allocatable :: values(:), matrix_values(:), all_values(:), all_matrix(:, :)
   integer :: i

   x = [(real(i, dp), i = 1, 10)]
   call check_close(moment(x), 5.5_dp, 1.0e-14_dp, 'first raw moment')
   call check_close(moment(x, 2, central=.true.), 8.25_dp, 1.0e-14_dp, 'second central moment')
   call check_close(moment(x, 1, central=.true., absolute=.true.), &
      2.5_dp, 1.0e-14_dp, 'mean absolute deviation')
   call check_close(skewness(x), 0.0_dp, 1.0e-14_dp, 'skewness')
   call check_close(kurtosis(x), 1.7757575757575756_dp, 1.0e-14_dp, 'kurtosis')
   call check_close(geary(x), 0.8703882797784892_dp, 1.0e-14_dp, 'geary')

   all_values = all_moments(x, 4)
   call check(size(all_values) == 5, 'all moments size')
   call check_close(all_values(1), 1.0_dp, 1.0e-14_dp, 'order zero')
   call check_close(all_values(2), 5.5_dp, 1.0e-14_dp, 'all moments order one')
   call check_close(all_values(3), 38.5_dp, 1.0e-14_dp, 'all moments order two')

   m(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   m(:, 2) = 2.0_dp * m(:, 1)
   matrix_values = moment(m, 2, central=.true.)
   call check_close(matrix_values(1), 2.0_dp, 1.0e-14_dp, 'matrix moment column one')
   call check_close(matrix_values(2), 8.0_dp, 1.0e-14_dp, 'matrix moment column two')
   values = skewness(m)
   call check(maxval(abs(values)) < 1.0e-14_dp, 'matrix skewness')
   all_matrix = all_moments(m, 3, central=.true.)
   call check_close(all_matrix(3, 1), 2.0_dp, 1.0e-14_dp, 'matrix all moments')

   xnan = [1.0_dp, 2.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 4.0_dp, 5.0_dp]
   call check_close(moment(xnan, na_rm=.true.), 3.0_dp, 1.0e-14_dp, 'na removal')

   print '(a)', 'test_moments: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      call check(abs(actual - expected) <= tolerance, label)
   end subroutine check_close

end program test_moments
