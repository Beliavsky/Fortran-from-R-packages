! SPDX-License-Identifier: GPL-3.0-or-later
module test_nonlinear_callbacks
   use cccp_kinds, only : dp
   implicit none
contains
   subroutine analytic_center(x,f,g,h,info)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::f,g(:),h(:,:)
      integer,intent(out)::info
      integer::j
      if(any(x<=0.0_dp))then
         info=1;f=huge(1.0_dp);g=0.0_dp;h=0.0_dp;return
      end if
      f=-sum(log(x));g=-1.0_dp/x;h=0.0_dp
      do j=1,size(x);h(j,j)=1.0_dp/(x(j)*x(j));end do
      info=0
   end subroutine analytic_center

   subroutine reciprocal_constraint(x,f,g,h,info)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::f(:),g(:,:),h(:,:,:)
      integer,intent(out)::info
      integer::j
      if(any(x<=0.0_dp))then
         info=1;f=huge(1.0_dp);g=0.0_dp;h=0.0_dp;return
      end if
      f(1)=sum(1.0_dp/x)-1.0_dp;g(1,:)=-1.0_dp/(x*x);h=0.0_dp
      do j=1,size(x);h(j,j,1)=2.0_dp/(x(j)**3);end do
      info=0
   end subroutine reciprocal_constraint
end module test_nonlinear_callbacks

program test_nonlinear
   use cccp
   use test_nonlinear_callbacks
   implicit none
   real(dp)::x0(3),a(1,3),b(1),q(2),x1(2),gmat(2,2),hv(2)
   type(cone_constraint)::c(1)
   type(cccp_solution)::sol
   integer::fails,i
   fails=0
   x0=0.25_dp;a=reshape([1.0_dp,1.0_dp,2.0_dp],[1,3]);b=1.0_dp
   call dcp(x0,analytic_center,a,b,control=ctrl(maxiters=150),sol=sol)
   call check(trim(sol%status)=='optimal','DCP status',fails)
   call check(abs(dot_product(a(1,:),sol%x)-1.0_dp)<1e-8_dp,'DCP equality',fails)
   call check(maxval(abs(sol%x-[1.0_dp/3.0_dp,1.0_dp/3.0_dp,1.0_dp/6.0_dp]))<2e-4_dp, &
      'DCP solution',fails)

   q=[1.0_dp,1.0_dp];x1=[1.5_dp,1.5_dp];gmat=0.0_dp
   do i=1,2;gmat(i,i)=-1.0_dp;end do
   hv=0.0_dp;c(1)=nnoc(gmat,hv)
   call dnl(q,x1,reciprocal_constraint,1,cones=c,control=ctrl(maxiters=150),sol=sol)
   call check(trim(sol%status)=='optimal','DNL status',fails)
   call check(maxval(abs(sol%x-[2.0_dp,2.0_dp]))<2e-3_dp,'DNL solution',fails)

   if(fails>0)error stop 1
   print '(a)','test_nonlinear: PASS'
contains
   subroutine check(ok,name,fails)
      logical,intent(in)::ok
      character(len=*),intent(in)::name
      integer,intent(inout)::fails
      if(.not.ok)then;print '(a,a)','FAIL: ',name;fails=fails+1;end if
   end subroutine check
end program test_nonlinear
