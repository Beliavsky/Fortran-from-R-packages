! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
module poilog_mle
   use poilog_kinds, only : dp
   use poilog_math, only : is_finite_dp
   use poilog_distribution, only : dpoilog, dbipoilog
   use poilog_rng, only : rpoilog, rbipoilog
   use poilog_optimize, only : optim_result, minimize
   implicit none
   private
   public :: poilog_fit, bipoilog_fit, poilog_mle_fit, bipoilog_mle_fit

   type :: poilog_fit
      real(dp) :: mu=0.0_dp, sig=0.0_dp, p=0.0_dp, loglik=-huge(1.0_dp)
      integer :: convergence=1, iterations=0
      real(dp), allocatable :: boot(:,:)
      real(dp) :: gof=-1.0_dp
   end type poilog_fit

   type :: bipoilog_fit
      real(dp) :: mu1=0.0_dp,mu2=0.0_dp,sig1=0.0_dp,sig2=0.0_dp,rho=0.0_dp
      real(dp) :: p(2)=0.0_dp, loglik=-huge(1.0_dp)
      integer :: convergence=1, iterations=0
      real(dp), allocatable :: boot(:,:)
      real(dp) :: gof=-1.0_dp
   end type bipoilog_fit

contains

   function poilog_mle_fit(n,start_mu,start_sig,nboot,ztrunc,method,maxit) result(fit)
      integer, intent(in) :: n(:)
      real(dp), intent(in), optional :: start_mu,start_sig
      integer, intent(in), optional :: nboot,maxit
      logical, intent(in), optional :: ztrunc
      character(len=*), intent(in), optional :: method
      type(poilog_fit) :: fit
      integer, allocatable :: un(:),nr(:),sim(:),bu(:),bnr(:)
      real(dp) :: x0(2),smu,ssig
      integer :: nb,mi,b
      logical :: zt
      character(len=32) :: meth
      type(optim_result) :: opt,bopt

      if(size(n)==0 .or. any(n<0)) return
      smu=1.0_dp; if(present(start_mu)) smu=start_mu
      ssig=2.0_dp; if(present(start_sig)) ssig=start_sig
      if(ssig<=0.0_dp) return
      nb=0; if(present(nboot)) nb=max(0,nboot)
      zt=.true.; if(present(ztrunc)) zt=ztrunc
      mi=1000; if(present(maxit)) mi=maxit
      meth='BFGS'; if(present(method)) meth=method
      call compress1(n,un,nr)
      x0=[smu,log(ssig)]
      opt=minimize(obj,x0,meth,mi,1.0e-6_dp)
      fit%mu=opt%par(1); fit%sig=exp(clamp_logsig(opt%par(2)))
      fit%p=1.0_dp-dpoilog(0,fit%mu,fit%sig)
      fit%loglik=-opt%value; fit%convergence=opt%convergence; fit%iterations=opt%iterations

      if(nb>0 .and. fit%convergence==0) then
         allocate(fit%boot(nb,3))
         do b=1,nb
            sim=rpoilog(size(n),fit%mu,fit%sig,cond_s=.true.,keep0=.not.zt)
            call compress1(sim,bu,bnr)
            bopt=minimize(bobj,opt%par,meth,mi,1.0e-6_dp)
            fit%boot(b,:)=[bopt%par(1),exp(clamp_logsig(bopt%par(2))),-bopt%value]
         end do
         fit%gof=real(1+count(fit%boot(:,3)<fit%loglik),dp)/real(nb,dp)
      end if

   contains
      function obj(z) result(v)
         real(dp), intent(in) :: z(:)
         real(dp) :: v
         v=nll1(z,un,nr,zt)
      end function obj
      function bobj(z) result(v)
         real(dp), intent(in) :: z(:)
         real(dp) :: v
         v=nll1(z,bu,bnr,zt)
      end function bobj
   end function poilog_mle_fit

   function bipoilog_mle_fit(n1,n2,start,nboot,ztrunc,method,maxit) result(fit)
      integer, intent(in) :: n1(:),n2(:)
      real(dp), intent(in), optional :: start(5)
      integer, intent(in), optional :: nboot,maxit
      logical, intent(in), optional :: ztrunc
      character(len=*), intent(in), optional :: method
      type(bipoilog_fit) :: fit
      integer, allocatable :: u1(:),u2(:),nr(:),sim(:,:),bu1(:),bu2(:),bnr(:)
      real(dp) :: x0(5),s(5)
      integer :: nb,mi,b
      logical :: zt
      character(len=32) :: meth
      type(optim_result) :: opt,bopt

      if(size(n1)==0 .or. size(n1)/=size(n2) .or. any(n1<0) .or. any(n2<0)) return
      s=[1.0_dp,1.0_dp,2.0_dp,2.0_dp,0.5_dp]; if(present(start)) s=start
      if(s(3)<=0.0_dp .or. s(4)<=0.0_dp .or. abs(s(5))>1.0_dp) return
      nb=0; if(present(nboot)) nb=max(0,nboot)
      zt=.true.; if(present(ztrunc)) zt=ztrunc
      mi=1000; if(present(maxit)) mi=maxit
      meth='BFGS'; if(present(method)) meth=method
      x0=[s(1),s(2),log(s(3)),log(s(4)),rho_transform(s(5))]
      call compress2(n1,n2,u1,u2,nr)
      opt=minimize(obj,x0,meth,mi,1.0e-5_dp)
      fit%mu1=opt%par(1); fit%mu2=opt%par(2)
      fit%sig1=exp(clamp_logsig(opt%par(3))); fit%sig2=exp(clamp_logsig(opt%par(4)))
      fit%rho=inv_rho_transform(opt%par(5))
      fit%p=[1.0_dp-dpoilog(0,fit%mu1,fit%sig1),1.0_dp-dpoilog(0,fit%mu2,fit%sig2)]
      fit%loglik=-opt%value; fit%convergence=opt%convergence; fit%iterations=opt%iterations
      if(nb>0 .and. fit%convergence==0) then
         allocate(fit%boot(nb,6))
         do b=1,nb
            sim=rbipoilog(size(n1),fit%mu1,fit%mu2,fit%sig1,fit%sig2,fit%rho,cond_s=.true.,keep0=.not.zt)
            call compress2(sim(:,1),sim(:,2),bu1,bu2,bnr)
            bopt=minimize(bobj,opt%par,meth,mi,1.0e-5_dp)
            fit%boot(b,:)=[bopt%par(1),bopt%par(2),exp(clamp_logsig(bopt%par(3))), &
               exp(clamp_logsig(bopt%par(4))),inv_rho_transform(bopt%par(5)),-bopt%value]
         end do
         fit%gof=real(1+count(fit%boot(:,6)<fit%loglik),dp)/real(nb,dp)
      end if
   contains
      function obj(z) result(v)
         real(dp), intent(in) :: z(:)
         real(dp) :: v
         v=nll2(z,u1,u2,nr,zt)
      end function obj
      function bobj(z) result(v)
         real(dp), intent(in) :: z(:)
         real(dp) :: v
         v=nll2(z,bu1,bu2,bnr,zt)
      end function bobj
   end function bipoilog_mle_fit

   function nll1(z,u,nr,zt) result(v)
      real(dp), intent(in) :: z(:)
      integer, intent(in) :: u(:),nr(:)
      logical, intent(in) :: zt
      real(dp) :: v,sig,bzero,p
      integer :: i
      sig=exp(clamp_logsig(z(2))); bzero=0.0_dp
      if(zt) then
         p=1.0_dp-dpoilog(0,z(1),sig)
         if(p<=0.0_dp) then; v=huge(1.0_dp); return; end if
         bzero=log(p)
      end if
      v=0.0_dp
      do i=1,size(u)
         p=dpoilog(u(i),z(1),sig)
         if(p<=0.0_dp .or. .not.is_finite_dp(p)) then; v=huge(1.0_dp); return; end if
         v=v-real(nr(i),dp)*(log(p)-bzero)
      end do
   end function nll1

   function nll2(z,u1,u2,nr,zt) result(v)
      real(dp), intent(in) :: z(:)
      integer, intent(in) :: u1(:),u2(:),nr(:)
      logical, intent(in) :: zt
      real(dp) :: v,sig1,sig2,rho,bzero,p
      integer :: i
      sig1=exp(clamp_logsig(z(3))); sig2=exp(clamp_logsig(z(4))); rho=inv_rho_transform(z(5))
      bzero=0.0_dp
      if(zt) then
         p=1.0_dp-dbipoilog(0,0,z(1),z(2),sig1,sig2,rho)
         if(p<=0.0_dp) then; v=huge(1.0_dp); return; end if
         bzero=log(p)
      end if
      v=0.0_dp
      do i=1,size(u1)
         p=dbipoilog(u1(i),u2(i),z(1),z(2),sig1,sig2,rho)
         if(p<=0.0_dp .or. .not.is_finite_dp(p)) then; v=huge(1.0_dp); return; end if
         v=v-real(nr(i),dp)*(log(p)-bzero)
      end do
   end function nll2

   subroutine compress1(x,u,nr)
      integer, intent(in) :: x(:)
      integer, allocatable, intent(out) :: u(:),nr(:)
      integer, allocatable :: tu(:),tn(:)
      integer :: i,j,m
      allocate(tu(size(x)),tn(size(x))); m=0; tn=0
      do i=1,size(x)
         j=0
         if(m>0) then
            do j=1,m; if(tu(j)==x(i)) exit; end do
            if(j>m) j=0
         end if
         if(j==0) then; m=m+1; tu(m)=x(i); tn(m)=1; else; tn(j)=tn(j)+1; end if
      end do
      allocate(u(m),nr(m)); u=tu(:m); nr=tn(:m)
   end subroutine compress1

   subroutine compress2(x,y,u1,u2,nr)
      integer, intent(in) :: x(:),y(:)
      integer, allocatable, intent(out) :: u1(:),u2(:),nr(:)
      integer, allocatable :: tx(:),ty(:),tn(:)
      integer :: i,j,m
      allocate(tx(size(x)),ty(size(x)),tn(size(x))); m=0; tn=0
      do i=1,size(x)
         j=0
         if(m>0) then
            do j=1,m; if(tx(j)==x(i).and.ty(j)==y(i)) exit; end do
            if(j>m) j=0
         end if
         if(j==0) then; m=m+1; tx(m)=x(i);ty(m)=y(i);tn(m)=1; else; tn(j)=tn(j)+1; end if
      end do
      allocate(u1(m),u2(m),nr(m)); u1=tx(:m);u2=ty(:m);nr=tn(:m)
   end subroutine compress2

   pure real(dp) function clamp_logsig(x) result(y)
      real(dp), intent(in) :: x
      y=min(354.0_dp,max(-372.0_dp,x))
   end function clamp_logsig

   pure real(dp) function rho_transform(rho) result(z)
      real(dp), intent(in) :: rho
      real(dp) :: r
      r=min(0.9999_dp,max(-0.9999_dp,rho))
      z=log((1.0_dp+r)/(1.0_dp-r))
   end function rho_transform

   pure real(dp) function inv_rho_transform(z) result(rho)
      real(dp), intent(in) :: z
      if(z>40.0_dp) then; rho=0.9999_dp
      else if(z<-40.0_dp) then; rho=-0.9999_dp
      else; rho=(1.0_dp-exp(-z))/(1.0_dp+exp(-z)); rho=min(0.9999_dp,max(-0.9999_dp,rho)); end if
   end function inv_rho_transform

end module poilog_mle
