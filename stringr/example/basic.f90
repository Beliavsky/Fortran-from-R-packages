! SPDX-License-Identifier: MIT
! SPDX-FileComment: Basic example for the Fortran stringr translation.
program basic
   !! Demonstrates elemental transformations and vector joining.
   use stringr
   implicit none
   character(len=12) :: names(3)
   names = ['  ALPHA     ', 'Beta value  ', 'gamma-value ']
   names = str_to_snake(names)
   print '(a)', str_flatten(names, ', ')
end program basic
