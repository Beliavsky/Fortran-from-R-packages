program test_latent
   use mcmcpack
   implicit none
   type(mcmc_result) :: r
   integer :: resp(5,3), md(8,4), yord(8)
   real(dp) :: th(5), al(3), be(3), ab0(2), abp(2,2), eq(5), iq(5)
   real(dp) :: pth(3), pal(2), peq(3), piq(3)
   real(dp) :: xfa(8,2), lam(2,1), psi(2), leq(2,1), lineq(2,1), lpm(2,1), lpp(2,1), a0(2), b0(2)
   real(dp) :: xo(8,2), bo(2), go(4), bp(2), bprec(2,2)
   integer :: i

   call set_seed(54321)
   resp = reshape([1,1,0,0,1, 0,1,1,0,1, 1,0,1,1,0],[5,3])
   th=0.0_dp; al=0.0_dp; be=1.0_dp; ab0=[0.0_dp,1.0_dp]
   abp=0.0_dp; abp(1,1)=0.5_dp; abp(2,2)=0.5_dp
   eq=-999.0_dp; eq(1)=0.0_dp; iq=0.0_dp
   r=mcmc_irt1d(resp,th,al,be,0.0_dp,1.0_dp,ab0,abp,eq,iq,10,20,2,.true.,.true.)
   call check(r%status==0 .and. size(r%draws,1)==10 .and. size(r%draws,2)==11,'irt1d')

   md=reshape([1,1,2,1, 1,1,3,3, 1,2,3,2, 1,2,1,1, &
               2,1,2,2, 2,1,3,1, 2,2,3,3, 2,3,1,1],[8,4],order=[2,1])
   ! Explicit assignment avoids any ambiguity in reshape order.
   md(1,:)=[1,1,2,1]; md(2,:)=[1,1,3,3]; md(3,:)=[1,2,3,2]; md(4,:)=[1,2,1,1]
   md(5,:)=[2,1,2,2]; md(6,:)=[2,1,3,1]; md(7,:)=[2,2,3,3]; md(8,:)=[2,3,1,1]
   pth=[-0.5_dp,0.0_dp,0.5_dp]; pal=[1.0_dp,1.0_dp]; peq=-999.0_dp; peq(2)=0.0_dp; piq=0.0_dp
   r=mcmc_paircompare(md,pth,pal,peq,piq,1.0_dp,1.0_dp,.false.,10,20,2,.true.,.true.)
   call check(r%status==0 .and. size(r%draws,1)==10 .and. size(r%draws,2)==5,'paircompare')

   do i=1,8
      xfa(i,1)=0.7_dp*real(i-4,dp)+0.1_dp*real(mod(i,2),dp)
      xfa(i,2)=0.5_dp*real(i-4,dp)-0.1_dp*real(mod(i,3),dp)
   end do
   lam(:,1)=[0.7_dp,0.5_dp]; psi=[1.0_dp,1.0_dp]; leq=-999.0_dp; lineq=0.0_dp
   lpm=0.0_dp; lpp=1.0_dp; a0=2.0_dp; b0=2.0_dp
   r=mcmc_factanal(xfa,lam,psi,leq,lineq,lpm,lpp,a0,b0,10,20,2,.true.)
   call check(r%status==0 .and. size(r%draws,1)==10 .and. size(r%draws,2)==12,'factanal')
   call check(all(r%draws(:,3:4)>0.0_dp),'factanal psi')

   xo(:,1)=1.0_dp
   xo(:,2)=[-1.5_dp,-1.0_dp,-0.5_dp,0.0_dp,0.2_dp,0.5_dp,1.0_dp,1.5_dp]
   yord=[1,1,1,2,2,2,3,3]; bo=0.0_dp; go=[-huge(1.0_dp)/10.0_dp,0.0_dp,0.8_dp,huge(1.0_dp)/10.0_dp]
   bp=0.0_dp; bprec=0.0_dp; bprec(1,1)=0.1_dp; bprec(2,2)=0.1_dp
   r=mcmc_oprobit(yord,xo,bo,go,bp,bprec,0.15_dp,10,20,2)
   call check(r%status==0 .and. size(r%draws,1)==10 .and. size(r%draws,2)==3,'oprobit')

   print '(a)', 'test_latent: PASS'
contains
   subroutine check(ok,label)
      logical,intent(in)::ok
      character(*),intent(in)::label
      if(.not.ok)then
         print '(a)', 'FAIL: '//label
         error stop 1
      end if
   end subroutine check
end program test_latent
