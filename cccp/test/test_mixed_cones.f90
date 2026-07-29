! SPDX-License-Identifier: GPL-3.0-or-later
program test_mixed_cones
   use cccp
   implicit none
   real(dp)::q(3),g(2,3),h(2),fmat(3,3),gv(3),d(3),fs
   real(dp)::fl(3,3,3),f0(3,3)
   type(cone_constraint)::c(3)
   type(cccp_solution)::sol
   integer::fails
   fails=0
   q=[-6.0_dp,-4.0_dp,-5.0_dp]
   g=row_matrix([16.0_dp,-14.0_dp,5.0_dp,7.0_dp,2.0_dp,0.0_dp],2,3)
   h=[-3.0_dp,5.0_dp];c(1)=nnoc(g,h)
   fmat=row_matrix([8.0_dp,13.0_dp,-12.0_dp,-8.0_dp,18.0_dp,6.0_dp,1.0_dp,-3.0_dp,-17.0_dp],3,3)
   gv=[-2.0_dp,-14.0_dp,-13.0_dp];d=[-24.0_dp,-7.0_dp,15.0_dp];fs=12.0_dp;c(2)=socc(fmat,gv,d,fs)
   fl(:,:,1)=reshape([7.0_dp,-5.0_dp,1.0_dp,-5.0_dp,1.0_dp,-7.0_dp,1.0_dp,-7.0_dp,-4.0_dp],[3,3])
   fl(:,:,2)=reshape([3.0_dp,13.0_dp,-6.0_dp,13.0_dp,12.0_dp,-10.0_dp,-6.0_dp,-10.0_dp,-28.0_dp],[3,3])
   fl(:,:,3)=reshape([9.0_dp,6.0_dp,-6.0_dp,6.0_dp,-7.0_dp,-7.0_dp,-6.0_dp,-7.0_dp,-11.0_dp],[3,3])
   f0=reshape([68.0_dp,-30.0_dp,-19.0_dp,-30.0_dp,99.0_dp,23.0_dp,-19.0_dp,23.0_dp,10.0_dp],[3,3])
   c(3)=psdc(fl,f0)
   call solve_lp(q,cones=c,control=ctrl(maxiters=180),sol=sol)
   call check(trim(sol%status)=='optimal','mixed status',fails)
   call check(sol%state%pslack>-2e-5_dp,'mixed feasibility',fails)
   if(fails>0)error stop 1
   print '(a)','test_mixed_cones: PASS'
contains
   function row_matrix(v,m,n) result(a)
      real(dp),intent(in)::v(:);integer,intent(in)::m,n;real(dp)::a(m,n)
      a=transpose(reshape(v,[n,m]))
   end function row_matrix
   subroutine check(ok,name,fails)
      logical,intent(in)::ok;character(len=*),intent(in)::name;integer,intent(inout)::fails
      if(.not.ok)then;print '(a,a)','FAIL: ',name;fails=fails+1;end if
   end subroutine check
end program test_mixed_cones
