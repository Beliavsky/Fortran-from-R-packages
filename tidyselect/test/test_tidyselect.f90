! SPDX-License-Identifier: MIT
! SPDX-FileComment: Tests for the Fortran tidyselect translation.
program test_tidyselect
   use tidyselect
   implicit none
   character(len=12) :: names(4)
   integer, allocatable :: selected(:)
   names=['id          ','return_us   ','return_eu   ','volume      ']
   selected=starts_with(names,'return_')
   if(any(selected/=[2,3])) error stop 'starts_with failed'
   if(any(ends_with(names,'ume')/=[4])) error stop 'ends_with failed'
   if(any(all_of(names,['id    ','volume'])/=[1,4])) error stop 'all_of failed'
   if(any(last_col(names)/=[4])) error stop 'last_col failed'
   if(any(eval_select(names,[-2,-3])/=[1,4])) error stop 'exclusion failed'
   print '(a)','All tidyselect tests passed.'
end program test_tidyselect
