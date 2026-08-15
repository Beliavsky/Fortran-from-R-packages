! General grouped random effects for GAMLSS parameters.
! Arbitrary random intercept/slope design with optional correlated effects.
! The covariance is updated by an EM-like posterior second-moment step while
! the GAMLSS coefficients are re-fitted by the RS/CG engine.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_random_effects_v03
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar, GAMLSS_NO
   use gamlss_types
   use gamlss_core, only : fit_gamlss_model
   use gamlss_linalg, only : invert_matrix
   use nlme_lme, only : fit_lme
   use nlme_types, only : lme_result
   implicit none
   private
   public :: random_effects_result_t, fit_gamlss_random_effects

   type, public :: random_effects_result_t
      type(gamlss_result_t) :: model
      real(dp), allocatable :: effects(:,:)
      real(dp), allocatable :: covariance(:,:), precision(:,:)
      integer, allocatable :: levels(:)
      integer :: parameter = 1
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
   end type random_effects_result_t
contains

   subroutine fit_gamlss_random_effects(y,x_fixed,z_random,group,family,result,parameter, &
      x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma,offset_nu,offset_tau,control, &
      correlated,use_nlme_start,max_outer,tol_cov)
      real(dp), intent(in) :: y(:),x_fixed(:,:),z_random(:,:)
      integer, intent(in) :: group(:),family
      type(random_effects_result_t), intent(out) :: result
      integer, intent(in), optional :: parameter,max_outer
      real(dp), intent(in), optional :: x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp), intent(in), optional :: offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      type(gamlss_control_t), intent(in), optional :: control
      logical, intent(in), optional :: correlated,use_nlme_start
      real(dp), intent(in), optional :: tol_cov

      real(dp), allocatable :: xm(:,:),xs(:,:),xn(:,:),xt(:,:),zbig(:,:)
      real(dp), allocatable :: pm(:,:),ps(:,:),pn(:,:),pt(:,:),sigma(:,:),sigma_new(:,:),prec(:,:)
      real(dp), allocatable :: start(:),effects(:,:),postcov(:,:),tmp(:,:)
      type(gamlss_control_t) :: ctl
      type(lme_result) :: lr
      integer :: n,np,target,q,ng,pfix,it,nouter,istat,g,a,b
      real(dp) :: crit,delta,eps
      logical :: corr,nlstart

      n=size(y);np=family_npar(family);target=1;if(present(parameter))target=parameter
      q=size(z_random,2)
      if(n<=0.or.size(x_fixed,1)/=n.or.size(z_random,1)/=n.or.size(group)/=n.or.q<1)then
         result%status=1;return
      end if
      if(target<1.or.target>np)then;result%status=2;return;end if
      call grouped_random_design(group,z_random,zbig,result%levels,istat)
      if(istat/=0)then;result%status=3;return;end if
      ng=size(result%levels)
      call default_design(n,x_fixed,x_sigma,x_nu,x_tau,xm,xs,xn,xt)
      select case(target)
      case(1);pfix=size(xm,2);call append_design(xm,zbig);call zero_penalty(xm,pm)
      case(2);pfix=size(xs,2);call append_design(xs,zbig);call zero_penalty(xs,ps)
      case(3);pfix=size(xn,2);call append_design(xn,zbig);call zero_penalty(xn,pn)
      case(4);pfix=size(xt,2);call append_design(xt,zbig);call zero_penalty(xt,pt)
      end select
      if(target/=1)call zero_penalty(xm,pm)
      if(target/=2)call zero_penalty(xs,ps)
      if(target/=3)call zero_penalty(xn,pn)
      if(target/=4)call zero_penalty(xt,pt)

      corr=.true.;if(present(correlated))corr=correlated
      nlstart=.true.;if(present(use_nlme_start))nlstart=use_nlme_start
      allocate(sigma(q,q));sigma=0.0_dp
      do a=1,q;sigma(a,a)=1.0_dp;end do
      if(nlstart.and.family==GAMLSS_NO.and.target==1)then
         call fit_lme(y,x_fixed,z_random,group,lr)
         if(lr%status==0.and.allocated(lr%random_covariance))then
            if(size(lr%random_covariance,1)==q.and.size(lr%random_covariance,2)==q)then
               sigma=lr%random_covariance
            end if
         end if
      end if
      if(.not.corr)call diagonalize(sigma)
      eps=1.0e-8_dp;call stabilize_covariance(sigma,eps)
      call invert_matrix(sigma,prec,istat)
      if(istat/=0)then;result%status=4;return;end if

      ctl=gamlss_control_t();if(present(control))ctl=control
      ! The full precision matrix is updated explicitly; do not also estimate a scalar lambda.
      select case(target)
      case(1);ctl%estimate_lambda_mu=.false.
      case(2);ctl%estimate_lambda_sigma=.false.
      case(3);ctl%estimate_lambda_nu=.false.
      case(4);ctl%estimate_lambda_tau=.false.
      end select
      nouter=20;if(present(max_outer))nouter=max(1,max_outer)
      crit=1.0e-5_dp;if(present(tol_cov))crit=max(0.0_dp,tol_cov)

      do it=1,nouter
         call install_random_penalty(target,pfix,ng,q,prec,pm,ps,pn,pt)
         if(it==1)then
            call fit_gamlss_model(y,xm,family,result%model,x_sigma=xs,x_nu=xn,x_tau=xt,weights=weights, &
               offset_mu=offset_mu,offset_sigma=offset_sigma,offset_nu=offset_nu,offset_tau=offset_tau, &
               penalty_mu=pm,penalty_sigma=ps,penalty_nu=pn,penalty_tau=pt, &
               lambda_mu=merge(1.0_dp,0.0_dp,target==1),lambda_sigma=merge(1.0_dp,0.0_dp,target==2), &
               lambda_nu=merge(1.0_dp,0.0_dp,target==3),lambda_tau=merge(1.0_dp,0.0_dp,target==4),control=ctl)
         else
            call flatten_start(result%model,np,start)
            call fit_gamlss_model(y,xm,family,result%model,x_sigma=xs,x_nu=xn,x_tau=xt,weights=weights, &
               offset_mu=offset_mu,offset_sigma=offset_sigma,offset_nu=offset_nu,offset_tau=offset_tau, &
               penalty_mu=pm,penalty_sigma=ps,penalty_nu=pn,penalty_tau=pt, &
               lambda_mu=merge(1.0_dp,0.0_dp,target==1),lambda_sigma=merge(1.0_dp,0.0_dp,target==2), &
               lambda_nu=merge(1.0_dp,0.0_dp,target==3),lambda_tau=merge(1.0_dp,0.0_dp,target==4), &
               start=start,control=ctl)
         end if
         if(result%model%status/=0)then;result%status=10+result%model%status;return;end if
         call extract_random_block(result%model,target,pfix,ng,q,effects,postcov,istat)
         if(istat/=0)then;result%status=5;return;end if
         allocate(sigma_new(q,q));sigma_new=0.0_dp
         do g=1,ng
            do a=1,q
               do b=1,q
                  sigma_new(a,b)=sigma_new(a,b)+effects(g,a)*effects(g,b)+ &
                     postcov((g-1)*q+a,(g-1)*q+b)
               end do
            end do
         end do
         sigma_new=sigma_new/real(ng,dp)
         if(.not.corr)call diagonalize(sigma_new)
         call stabilize_covariance(sigma_new,eps)
         delta=maxval(abs(sigma_new-sigma)/(1.0_dp+abs(sigma)))
         tmp=0.5_dp*sigma+0.5_dp*sigma_new
         call move_alloc(tmp,sigma)
         deallocate(sigma_new)
         call invert_matrix(sigma,prec,istat)
         if(istat/=0)then;result%status=6;return;end if
         result%iterations=it
         if(delta<crit)then;result%converged=.true.;exit;end if
      end do
      result%effects=effects;result%covariance=sigma;result%precision=prec
      result%parameter=target;result%status=0
   contains
      subroutine default_design(nn,xf,xs_in,xn_in,xt_in,a1,a2,a3,a4)
         integer,intent(in)::nn
         real(dp),intent(in)::xf(:,:)
         real(dp),intent(in),optional::xs_in(:,:),xn_in(:,:),xt_in(:,:)
         real(dp),allocatable,intent(out)::a1(:,:),a2(:,:),a3(:,:),a4(:,:)
         a1=xf
         if(present(xs_in))then;a2=xs_in;else;allocate(a2(nn,1));a2=1.0_dp;end if
         if(present(xn_in))then;a3=xn_in;else;allocate(a3(nn,1));a3=1.0_dp;end if
         if(present(xt_in))then;a4=xt_in;else;allocate(a4(nn,1));a4=1.0_dp;end if
      end subroutine default_design
   end subroutine fit_gamlss_random_effects

   subroutine grouped_random_design(group,z,zbig,levels,status)
      integer,intent(in)::group(:)
      real(dp),intent(in)::z(:,:)
      real(dp),allocatable,intent(out)::zbig(:,:)
      integer,allocatable,intent(out)::levels(:)
      integer,intent(out)::status
      integer,allocatable::tmp(:)
      integer::i,g,q,ng,idx
      status=0;q=size(z,2)
      if(size(group)/=size(z,1).or.size(group)==0.or.q<1)then
         allocate(zbig(0,0),levels(0));status=1;return
      end if
      allocate(tmp(size(group)));ng=0
      do i=1,size(group)
         if(ng==0)then
            ng=1;tmp(1)=group(i)
         else if(.not.any(tmp(1:ng)==group(i)))then
            ng=ng+1;tmp(ng)=group(i)
         end if
      end do
      allocate(levels(ng));levels=tmp(1:ng)
      allocate(zbig(size(group),ng*q));zbig=0.0_dp
      do i=1,size(group)
         g=0
         do idx=1,ng
            if(group(i)==levels(idx))then;g=idx;exit;end if
         end do
         zbig(i,(g-1)*q+1:g*q)=z(i,:)
      end do
   end subroutine grouped_random_design

   subroutine append_design(x,z)
      real(dp),allocatable,intent(inout)::x(:,:)
      real(dp),intent(in)::z(:,:)
      real(dp),allocatable::tmp(:,:)
      integer::p
      p=size(x,2);allocate(tmp(size(x,1),p+size(z,2)))
      tmp(:,1:p)=x;tmp(:,p+1:)=z;call move_alloc(tmp,x)
   end subroutine append_design

   subroutine zero_penalty(x,p)
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable,intent(out)::p(:,:)
      allocate(p(size(x,2),size(x,2)));p=0.0_dp
   end subroutine zero_penalty

   subroutine install_random_penalty(target,pfix,ng,q,prec,pm,ps,pn,pt)
      integer,intent(in)::target,pfix,ng,q
      real(dp),intent(in)::prec(:,:)
      real(dp),intent(inout)::pm(:,:),ps(:,:),pn(:,:),pt(:,:)
      integer::g,lo,hi
      select case(target)
      case(1);pm=0.0_dp
      case(2);ps=0.0_dp
      case(3);pn=0.0_dp
      case(4);pt=0.0_dp
      end select
      do g=1,ng
         lo=pfix+(g-1)*q+1;hi=lo+q-1
         select case(target)
         case(1);pm(lo:hi,lo:hi)=prec
         case(2);ps(lo:hi,lo:hi)=prec
         case(3);pn(lo:hi,lo:hi)=prec
         case(4);pt(lo:hi,lo:hi)=prec
         end select
      end do
   end subroutine install_random_penalty

   subroutine extract_random_block(model,target,pfix,ng,q,effects,postcov,status)
      type(gamlss_result_t),intent(in)::model
      integer,intent(in)::target,pfix,ng,q
      real(dp),allocatable,intent(out)::effects(:,:),postcov(:,:)
      integer,intent(out)::status
      real(dp),allocatable::coef(:),cov(:,:)
      integer::g,lo,hi
      status=0
      select case(target)
      case(1);coef=model%mu%coefficients;cov=model%mu%covariance
      case(2);coef=model%sigma%coefficients;cov=model%sigma%covariance
      case(3);coef=model%nu%coefficients;cov=model%nu%covariance
      case(4);coef=model%tau%coefficients;cov=model%tau%covariance
      case default;status=1;return
      end select
      if(size(coef)<pfix+ng*q.or.size(cov,1)<pfix+ng*q)then;status=2;return;end if
      allocate(effects(ng,q),postcov(ng*q,ng*q));postcov=0.0_dp
      do g=1,ng
         lo=pfix+(g-1)*q+1;hi=lo+q-1
         effects(g,:)=coef(lo:hi)
         postcov((g-1)*q+1:g*q,(g-1)*q+1:g*q)=cov(lo:hi,lo:hi)
      end do
   end subroutine extract_random_block

   subroutine flatten_start(model,np,start)
      type(gamlss_result_t),intent(in)::model
      integer,intent(in)::np
      real(dp),allocatable,intent(out)::start(:)
      integer::n1,n2,n3,n4,pos
      n1=size(model%mu%coefficients);n2=0;n3=0;n4=0
      if(np>=2)n2=size(model%sigma%coefficients)
      if(np>=3)n3=size(model%nu%coefficients)
      if(np>=4)n4=size(model%tau%coefficients)
      allocate(start(n1+n2+n3+n4));pos=0
      start(1:n1)=model%mu%coefficients;pos=n1
      if(np>=2)then;start(pos+1:pos+n2)=model%sigma%coefficients;pos=pos+n2;end if
      if(np>=3)then;start(pos+1:pos+n3)=model%nu%coefficients;pos=pos+n3;end if
      if(np>=4)start(pos+1:pos+n4)=model%tau%coefficients
   end subroutine flatten_start

   subroutine diagonalize(a)
      real(dp),intent(inout)::a(:,:)
      integer::i,j
      do i=1,size(a,1);do j=1,size(a,2);if(i/=j)a(i,j)=0.0_dp;end do;end do
   end subroutine diagonalize

   subroutine stabilize_covariance(a,eps)
      real(dp),intent(inout)::a(:,:)
      real(dp),intent(in)::eps
      integer::i
      a=0.5_dp*(a+transpose(a))
      do i=1,size(a,1);a(i,i)=max(a(i,i),eps);end do
   end subroutine stabilize_covariance

end module gamlss_random_effects_v03
