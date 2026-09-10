! SPDX-License-Identifier: MIT
! SPDX-FileComment: Map-reduce example for the Fortran purrr translation.
program map_reduce
   !! Maps a square callback and reduces with addition.
   use purrr
   implicit none
   integer, allocatable :: squares(:)
   squares = map([1, 2, 3, 4], square)
   print '(a,*(1x,i0))', 'Squares:', squares
   print '(a,1x,i0)', 'Sum:', reduce(squares, add)
contains
   pure integer function square(value)
      !! Squares one integer.
      integer, intent(in) :: value !! Value to square.
      square = value * value
   end function square
   pure integer function add(left, right)
      !! Adds two integers.
      integer, intent(in) :: left  !! Left addend.
      integer, intent(in) :: right !! Right addend.
      add = left + right
   end function add
end program map_reduce
