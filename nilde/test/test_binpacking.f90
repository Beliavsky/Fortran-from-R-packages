! SPDX-License-Identifier: GPL-2.0-or-later
program test_binpacking
   use nilde
   implicit none
   type(bin_packing_result_t) :: r
   integer(i8) :: w(7)
   integer :: j,b
   w=[70_i8,60_i8,50_i8,40_i8,30_i8,20_i8,10_i8]
   r=bin_packing(w,100_i8)
   call check(r%min_bins==3,'minimum bins')
   call check(r%nsol==9,'optimal packing count')
   do j=1,r%nsol
      do b=1,r%min_bins
         call check(sum(w,mask=r%assignment(:,j)==b)<=100_i8,'bin capacity')
      end do
      call check(r%total_ineff(j)==20_i8,'total inefficiency')
   end do
   print *, 'test_binpacking: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(len=*),intent(in)::msg
      if(.not.ok) then; print *,'FAIL: ',trim(msg); error stop 1; end if
   end subroutine
end program
