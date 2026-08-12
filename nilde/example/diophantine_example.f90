! SPDX-License-Identifier: GPL-2.0-or-later
program diophantine_example
   use nilde
   implicit none
   type(integer_solutions_t) :: r
   integer :: j

   r = nlde([3_i8,2_i8,5_i8,16_i8], 18_i8, m=6, at_most=.false.)
   print '(a,i0)', 'number of solutions = ', r%nsol
   do j = 1, r%nsol
      print '(*(i0,1x))', r%x(:,j)
   end do
end program diophantine_example
