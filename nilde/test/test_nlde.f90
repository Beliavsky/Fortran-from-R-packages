! SPDX-License-Identifier: GPL-2.0-or-later
program test_nlde
   use nilde
   implicit none
   type(integer_solutions_t) :: r
   integer(i8) :: a4(4), a2(2), a9(9)
   integer :: j

   a4=[3_i8,2_i8,5_i8,16_i8]
   r=nlde(a4,18_i8)
   call check(r%nsol==10,'nlde count')
   do j=1,r%nsol
      call check(sum(a4*r%x(:,j))==18_i8,'nlde equation')
   end do

   r=nlde(a4,18_i8,m=6,at_most=.false.)
   call check(r%nsol==3,'exact M count')
   do j=1,r%nsol
      call check(sum(r%x(:,j))==6_i8,'exact M')
   end do

   r=nlde(a4,18_i8,m=6,option=1)
   call check(r%nsol==1,'binary nlde count')
   call check(all(r%x(:,1)==[0_i8,1_i8,0_i8,1_i8]),'binary nlde solution')

   a2=[15_i8,21_i8]
   r=nlde(a2,261_i8)
   call check(r%nsol==3,'two coefficient count')

   a9=[70_i8,60_i8,50_i8,33_i8,33_i8,33_i8,11_i8,7_i8,3_i8]
   r=nlde(a9,100_i8,option=2)
   call check(r%nsol==108,'binary inequality count')
   do j=1,r%nsol
      call check(sum(a9*r%x(:,j))<=100_i8,'binary inequality')
      call check(all(r%x(:,j)==0_i8 .or. r%x(:,j)==1_i8),'binary inequality values')
   end do
   print *, 'test_nlde: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(len=*),intent(in)::msg
      if(.not.ok) then; print *,'FAIL: ',trim(msg); error stop 1; end if
   end subroutine
end program
