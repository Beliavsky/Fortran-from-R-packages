! SPDX-License-Identifier: GPL-2.0-or-later
program test_tsp
   use nilde
   implicit none
   type(tsp_result_t) :: r
   integer(i8) :: c(4,4)
   integer :: j,k
   integer(i8) :: val
   c=reshape([0_i8,1_i8,10_i8,10_i8, &
              1_i8,0_i8,10_i8,10_i8, &
              10_i8,10_i8,0_i8,1_i8, &
              10_i8,10_i8,1_i8,0_i8],[4,4])
   ! Matrix is symmetric, so reshape orientation is immaterial here.
   call check(assignment_lower_bound(c+diag_penalty())==4_i8,'assignment lower bound')
   r=tsp_solver(c)
   call check(r%initial_lower_bound==4_i8,'tsp initial lower')
   call check(r%tour_length==22_i8,'tsp optimum')
   call check(r%iterations==19,'tsp lower-bound iterations')
   call check(r%ntours==4,'tsp optimal tours')
   do j=1,r%ntours
      val=0_i8
      do k=1,4
         val=val+c(r%tours(k,j),r%tours(merge(k+1,1,k<4),j))
      end do
      call check(val==22_i8,'tour cost')
   end do
   print *, 'test_tsp: PASS'
contains
   function diag_penalty() result(p)
      integer(i8) :: p(4,4)
      integer :: i
      p=0_i8
      do i=1,4; p(i,i)=1000000_i8; end do
   end function
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(len=*),intent(in)::msg
      if(.not.ok) then; print *,'FAIL: ',trim(msg); error stop 1; end if
   end subroutine
end program
