program test_parity_v03
   use directional, only : dp, desag, desag3, resag, distance_cor, circ_dcor, spher_dcor
   use directional, only : knn_reg, knnreg_tune, rspcauchy, rpkbd, mixspcauchy_mle, mixpkbd_mle
   use directional, only : dmixspcauchy, dmixpkbd, rmixspcauchy, rmixpkbd, rkent, kent_mle
   implicit none
   integer :: fails
   fails=0
   call test_esag(fails)
   call test_dcor(fails)
   call test_knnreg(fails)
   call test_mixtures(fails)
   call test_kent(fails)
   if(fails==0)then
      print '(a)', 'test_parity_v03: PASS'
   else
      print '(a,i0)', 'test_parity_v03: FAIL ',fails
      error stop 1
   end if
contains
   subroutine check(ok,fails)
      logical,intent(in)::ok;integer,intent(inout)::fails;if(.not.ok)fails=fails+1
   end subroutine
   subroutine test_esag(fails)
      integer,intent(inout)::fails
      real(dp)::mu3(3),gam3(2),y3(3,3),a(3),b(3),mu4(4),gam4(5),y4(30,4),d4(30),nr
      integer::i
      mu3=[1.0_dp,0.7_dp,-0.3_dp];gam3=0.0_dp
      y3=reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp],[3,3])
      a=desag(y3,mu3,gam3);b=desag3(y3,mu3,gam3)
      call check(maxval(abs(a-b))<1e-10_dp,fails)
      mu4=[0.7_dp,-0.2_dp,0.4_dp,0.8_dp];gam4=[0.2_dp,-0.1_dp,0.1_dp,0.05_dp,-0.12_dp]
      y4=resag(30,mu4,gam4);d4=desag(y4,mu4,gam4)
      do i=1,30;nr=sqrt(sum(y4(i,:)**2));call check(abs(nr-1.0_dp)<1e-10_dp,fails);end do
      call check(all(d4>0.0_dp).and.all(d4<huge(1.0_dp)),fails)
   end subroutine
   subroutine test_dcor(fails)
      integer,intent(inout)::fails
      real(dp)::x(5,2),th(5),r
      x=reshape([0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp, 1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp],[5,2])
      r=distance_cor(x,x);call check(abs(r-1.0_dp)<1e-10_dp,fails)
      th=[0.0_dp,30.0_dp,80.0_dp,160.0_dp,250.0_dp];r=circ_dcor(th,th);call check(abs(r-1.0_dp)<1e-10_dp,fails)
      r=spher_dcor(x,x);call check(abs(r-1.0_dp)<1e-10_dp,fails)
   end subroutine
   subroutine test_knnreg(fails)
      integer,intent(inout)::fails
      real(dp)::x(8,1),y(8,1),xn(2,1),est(2,1),crit(4)
      integer::i,best
      do i=1,8;x(i,1)=real(i,dp);y(i,1)=2.0_dp*x(i,1);end do
      xn(:,1)=[2.1_dp,6.9_dp];call knn_reg(xn,y,x,1,est)
      call check(abs(est(1,1)-4.0_dp)<1e-12_dp.and.abs(est(2,1)-14.0_dp)<1e-12_dp,fails)
      call knnreg_tune(y,x,5,best,crit,nfolds=4);call check(best>=2.and.best<=5.and.all(crit>=0.0_dp),fails)
   end subroutine
   subroutine test_mixtures(fails)
      integer,intent(inout)::fails
      real(dp)::x1(50,3),x2(50,3),x(100,3),pr(2),rho(2),mu(2,3),ll,d(100),s,xr(20,3)
      integer::cl(100),idr(20),i
      x1=rspcauchy(50,[1.0_dp,0.0_dp,0.0_dp],0.75_dp);x2=rspcauchy(50,[-1.0_dp,0.0_dp,0.0_dp],0.75_dp)
      x(1:50,:)=x1;x(51:100,:)=x2
      call mixspcauchy_mle(x,2,pr,rho,mu,ll,cl,maxiters=30)
      call check(abs(sum(pr)-1.0_dp)<1e-8_dp.and.all(rho>=0.0_dp).and.all(rho<1.0_dp),fails)
      d=dmixspcauchy(x,pr,mu,rho);call check(all(d>0.0_dp),fails)
      call mixpkbd_mle(x,2,pr,rho,mu,ll,cl,maxiters=20)
      call check(abs(sum(pr)-1.0_dp)<1e-8_dp.and.all(rho>=0.0_dp).and.all(rho<1.0_dp),fails)
      d=dmixpkbd(x,pr,mu,rho);s=sum(d);call check(s>0.0_dp,fails)
      xr=rpkbd(20,[1.0_dp,0.0_dp,0.0_dp],0.5_dp);do i=1,20;call check(abs(sum(xr(i,:)**2)-1.0_dp)<1e-10_dp,fails);end do
      call rmixspcauchy(20,[0.4_dp,0.6_dp],reshape([1.0_dp,0.0_dp,0.0_dp,-1.0_dp,0.0_dp,0.0_dp],[2,3]),[0.4_dp,0.6_dp],xr,idr)
      call check(all(idr>=1).and.all(idr<=2),fails)
      call rmixpkbd(20,[0.5_dp,0.5_dp],reshape([1.0_dp,0.0_dp,0.0_dp,-1.0_dp,0.0_dp,0.0_dp],[2,3]),[0.3_dp,0.5_dp],xr,idr)
      call check(all(idr>=1).and.all(idr<=2),fails)
   end subroutine
   subroutine test_kent(fails)
      integer,intent(inout)::fails
      real(dp)::x(80,3),g(3,3),k,b,psi,lc,ll
      x=rkent(80,4.0_dp,[1.0_dp,0.0_dp,0.0_dp],0.7_dp)
      call kent_mle(x,g,k,b,psi,lc,ll)
      call check(k>0.0_dp.and.b>=0.0_dp.and.2.0_dp*b<k,fails)
      call check(abs(sum(g(:,1)*g(:,1))-1.0_dp)<1e-8_dp,fails)
      call check(ll>-huge(1.0_dp),fails)
   end subroutine
end program
