! SPDX-License-Identifier: GPL-2.0-or-later
module gend_example_functions
   use numderiv, only : dp
   implicit none
contains
   function model(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: value(:)
      value = [x(1), x(1) * x(2), x(2) ** 2]
   end function model
end module gend_example_functions

program gend_example
   use gend_example_functions, only : model
   use numderiv, only : dp, gend, gend_result
   implicit none
   type(gend_result) :: result
   integer :: i

   call gend(model, [2.0_dp, 3.0_dp], result)
   print '(a)', 'Columns: dx1, dx2, d2x1x1, d2x2x1, d2x2x2'
   do i = 1, size(result%dmat, 1)
      print '(*(1x,es16.8))', result%dmat(i, :)
   end do
end program gend_example
