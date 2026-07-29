! SPDX-License-Identifier: GPL-3.0-or-later
program test_cones
   use cccp
   implicit none
   real(dp) :: q(3), f1(2,3), g1(2), d1(3), fsc1
   real(dp) :: f2(3,3), g2(3), d2(3), fsc2
   type(cone_constraint) :: cones2(2), cones_psd(2)
   type(cccp_solution) :: sol
   real(dp) :: fp1(2,2,3), f01(2,2), fp2(3,3,3), f02(3,3)
   integer :: fails
   fails=0

   q=[-2.0_dp,1.0_dp,5.0_dp]
   f1=row_matrix([-13.0_dp,3.0_dp,5.0_dp,-12.0_dp,12.0_dp,-6.0_dp],2,3)
   g1=[-3.0_dp,-2.0_dp];d1=[-12.0_dp,-6.0_dp,5.0_dp];fsc1=-12.0_dp
   f2=row_matrix([-3.0_dp,6.0_dp,2.0_dp,1.0_dp,9.0_dp,2.0_dp,-1.0_dp,-19.0_dp,3.0_dp],3,3)
   g2=[0.0_dp,3.0_dp,-42.0_dp];d2=[-3.0_dp,6.0_dp,-10.0_dp];fsc2=27.0_dp
   cones2(1)=socc(f1,g1,d1,fsc1);cones2(2)=socc(f2,g2,d2,fsc2)
   call solve_lp(q,cones=cones2,control=ctrl(maxiters=150,feastol=1e-6_dp),sol=sol)
   call check(trim(sol%status)=='optimal','SOCP status',fails)
   call check(soc_feasible(cones2(1),sol%x,2e-5_dp).and.soc_feasible(cones2(2),sol%x,2e-5_dp), &
      'SOCP feasible',fails)

   q=[1.0_dp,-1.0_dp,1.0_dp]
   fp1(:,:,1)=reshape([-7.0_dp,-11.0_dp,-11.0_dp,3.0_dp],[2,2])
   fp1(:,:,2)=reshape([7.0_dp,-18.0_dp,-18.0_dp,8.0_dp],[2,2])
   fp1(:,:,3)=reshape([-2.0_dp,-8.0_dp,-8.0_dp,1.0_dp],[2,2])
   f01=reshape([33.0_dp,-9.0_dp,-9.0_dp,26.0_dp],[2,2])
   fp2(:,:,1)=reshape([-21.0_dp,-11.0_dp,0.0_dp,-11.0_dp,10.0_dp,8.0_dp,0.0_dp,8.0_dp,5.0_dp],[3,3])
   fp2(:,:,2)=reshape([0.0_dp,10.0_dp,16.0_dp,10.0_dp,-10.0_dp,-10.0_dp,16.0_dp,-10.0_dp,3.0_dp],[3,3])
   fp2(:,:,3)=reshape([-5.0_dp,2.0_dp,-17.0_dp,2.0_dp,-6.0_dp,8.0_dp,-17.0_dp,8.0_dp,6.0_dp],[3,3])
   f02=reshape([14.0_dp,9.0_dp,40.0_dp,9.0_dp,91.0_dp,10.0_dp,40.0_dp,10.0_dp,15.0_dp],[3,3])
   cones_psd(1)=psdc(fp1,f01);cones_psd(2)=psdc(fp2,f02)
   call solve_lp(q,cones=cones_psd,control=ctrl(maxiters=180),sol=sol)
   call check(trim(sol%status)=='optimal','SDP status',fails)
   call check(sol%state%pslack>-5e-5_dp,'SDP feasible',fails)

   if(fails>0)error stop 1
   print '(a)','test_cones: PASS'
contains
   function row_matrix(v,m,n) result(a)
      real(dp),intent(in)::v(:)
      integer,intent(in)::m,n
      real(dp)::a(m,n)
      a=transpose(reshape(v,[n,m]))
   end function row_matrix
   logical function soc_feasible(c,x,tol) result(ok)
      type(cone_constraint),intent(in)::c
      real(dp),intent(in)::x(:),tol
      real(dp)::s(size(c%h))
      s=c%h-matmul(c%g,x)
      ok=s(1)>=sqrt(dot_product(s(2:),s(2:)))-tol
   end function soc_feasible
   subroutine check(ok,name,fails)
      logical,intent(in)::ok
      character(len=*),intent(in)::name
      integer,intent(inout)::fails
      if(.not.ok)then;print '(a,a)','FAIL: ',name;fails=fails+1;end if
   end subroutine check
end program test_cones
