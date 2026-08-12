! SPDX-License-Identifier: GPL-2.0-or-later
program test_partitions_subsetsum
   use nilde
   implicit none
   type(integer_solutions_t) :: r
   integer(i8) :: a10(10), a5(5), b5(5)
   integer :: j

   r=get_partitions(8_i8,6,.true.)
   call check(r%nsol==20,'partitions at most count')
   do j=1,r%nsol
      call check(sum(r%x(:,j))==8_i8,'partition sum')
      call check(all(r%x(2:,j)>=r%x(:5,j)),'partition order')
   end do
   r=get_partitions(8_i8,6,.false.)
   call check(r%nsol==2,'partitions exactly count')
   call check(all(r%x>=1_i8),'exact partitions positive')

   a10=[41_i8,34_i8,21_i8,20_i8,8_i8,7_i8,7_i8,4_i8,3_i8,3_i8]
   r=get_subsetsum(a10,50_i8,m=10,problem='subsetsum01')
   call check(r%nsol==2,'subset01 count')
   do j=1,r%nsol
      call check(sum(a10*r%x(:,j))==50_i8,'subset01 sum')
   end do

   a5=[30_i8,29_i8,32_i8,31_i8,33_i8]
   b5=[1_i8,2_i8,1_i8,3_i8,4_i8]
   r=get_subsetsum(a5,91_i8,m=5,problem='bsubsetsum',bounds=b5)
   call check(r%nsol==3,'bounded subset count')
   do j=1,r%nsol
      call check(sum(a5*r%x(:,j))==91_i8,'bounded subset sum')
      call check(all(r%x(:,j)<=b5),'bounded subset bounds')
   end do
   print *, 'test_partitions_subsetsum: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(len=*),intent(in)::msg
      if(.not.ok) then; print *,'FAIL: ',trim(msg); error stop 1; end if
   end subroutine
end program
