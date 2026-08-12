! SPDX-License-Identifier: GPL-2.0-or-later
program test_knapsack
   use nilde
   implicit none
   type(knapsack_result_t) :: r
   real(dp) :: obj3(3)
   integer(i8) :: a3(3), b3(3)
   integer :: j

   obj3=[4.0_dp,3.0_dp,3.0_dp]; a3=[3_i8,2_i8,2_i8]; b3=[2_i8,2_i8,2_i8]
   r=get_knapsack(obj3,a3,4_i8,problem='bknap',bounds=b3)
   call check(abs(r%objective-6.0_dp)<1e-12_dp,'bknap objective')
   call check(r%nsol==3,'bknap tie count')
   do j=1,r%nsol
      call check(sum(a3*r%x(:,j))<=4_i8,'bknap capacity')
   end do

   r=get_knapsack(obj3,a3,4_i8,problem='knap01')
   call check(abs(r%objective-6.0_dp)<1e-12_dp,'knap01 objective')
   call check(r%nsol==1,'knap01 count')
   call check(all(r%x(:,1)==[0_i8,1_i8,1_i8]),'knap01 solution')

   r=get_knapsack([5.0_dp,6.0_dp],[2_i8,3_i8],6_i8,problem='uknap')
   call check(r%legacy_unbounded_all,'legacy unbounded flag')
   call check(r%nsol>1,'legacy unbounded enumerates feasible vectors')
   r=get_knapsack([5.0_dp,6.0_dp],[2_i8,3_i8],6_i8,problem='uknap',legacy_unbounded_all=.false.)
   call check(abs(r%objective-15.0_dp)<1e-12_dp,'corrected unbounded objective')
   call check(r%nsol==1 .and. all(r%x(:,1)==[3_i8,0_i8]),'corrected unbounded solution')
   print *, 'test_knapsack: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(len=*),intent(in)::msg
      if(.not.ok) then; print *,'FAIL: ',trim(msg); error stop 1; end if
   end subroutine
end program
