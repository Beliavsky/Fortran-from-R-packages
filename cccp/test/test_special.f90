! SPDX-License-Identifier: GPL-3.0-or-later
program test_special
   use cccp
   implicit none
   real(dp)::p(4,2),q(4),cov(3,3),x0(3),mrc(3),f0(3,2),g0(3),f1(1,2),g1(1)
   type(cccp_solution)::sol
   type(gp_function)::gc(1)
   integer::fails
   fails=0
   p=reshape([1.0_dp,1.0_dp,1.0_dp,1.0_dp, 0.0_dp,1.0_dp,2.0_dp,3.0_dp],[4,2])
   q=[1.0_dp,2.0_dp,2.0_dp,4.0_dp]
   sol=l1(p,q)
   call check(trim(sol%status)=='optimal','L1 status',fails)
   call check(sum(abs(matmul(p,sol%x(1:2))-q))<1.01_dp,'L1 objective',fails)

   cov=reshape([0.04_dp,0.006_dp,0.004_dp,0.006_dp,0.09_dp,0.01_dp,0.004_dp,0.01_dp,0.16_dp],[3,3])
   x0=[0.4_dp,0.35_dp,0.25_dp];mrc=[1.0_dp,1.0_dp,1.0_dp]/3.0_dp
   sol=rp(x0,cov,mrc)
   call check(trim(sol%status)=='optimal','RP status',fails)
   call check(abs(sum(sol%x)-1.0_dp)<1e-10_dp.and.all(sol%x>0.0_dp),'RP weights',fails)
   call check(maxval(abs(sol%x*matmul(cov,sol%x)-sum(sol%x*matmul(cov,sol%x))/3.0_dp))<2e-5_dp, &
      'RP contributions',fails)

   f0=transpose(reshape([3.0_dp,-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,-3.0_dp],[2,3]))
   g0=log([0.44_dp,10.0_dp,0.592_dp])
   f1=reshape([-1.0_dp,3.0_dp],[1,2]);g1=log([8.62_dp])
   allocate(gc(1)%f(1,2),gc(1)%g(1));gc(1)%f=f1;gc(1)%g=g1
   sol=gp(f0,g0,gc)
   call check(trim(sol%status)=='optimal','GP status',fails)
   call check(all(sol%x>0.0_dp),'GP positive',fails)
   call check(maxval(abs(sol%x-[1.2866774442538667_dp,0.5304618317646771_dp]))<2e-5_dp, &
      'GP reference',fails)

   if(fails>0)error stop 1
   print '(a)','test_special: PASS'
contains
   subroutine check(ok,name,fails)
      logical,intent(in)::ok
      character(len=*),intent(in)::name
      integer,intent(inout)::fails
      if(.not.ok)then;print '(a,a)','FAIL: ',name;fails=fails+1;end if
   end subroutine check
end program test_special
