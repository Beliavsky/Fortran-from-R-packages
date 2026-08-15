! Simultaneous arbitrary grouped random-effect blocks on multiple GAMLSS parameters.
! Each active parameter has its own q-dimensional within-group covariance matrix.
! Cross-parameter random-effect covariances are intentionally not imposed here.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_multi_random_v05
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar
   use gamlss_types
   use gamlss_core, only : fit_gamlss_model
   use gamlss_linalg, only : invert_matrix
   implicit none
   private
   public :: multi_random_effects_result_t,fit_gamlss_multi_random_effects

   type,public :: multi_random_effects_result_t
      type(gamlss_result_t) :: model
      real(dp),allocatable :: effects(:,:,:)
      real(dp),allocatable :: covariance(:,:,:)
      real(dp),allocatable :: precision(:,:,:)
      integer,allocatable :: levels(:)
      logical :: active(4)=.false.
      logical :: correlated_within(4)=.true.
      integer :: iterations=0
      integer :: status=0
      logical :: converged=.false.
   end type multi_random_effects_result_t
contains

   subroutine fit_gamlss_multi_random_effects(y,x_mu,z_random,group,family,result,active_parameters, &
      correlated_within,x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma,offset_nu,offset_tau, &
      control,max_outer,tol_cov)
      real(dp),intent(in)::y(:),x_mu(:,:),z_random(:,:,:)
      integer,intent(in)::group(:),family
      type(multi_random_effects_result_t),intent(out)::result
      logical,intent(in),optional::active_parameters(4),correlated_within(4)
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      type(gamlss_control_t),intent(in),optional::control
      integer,intent(in),optional::max_outer
      real(dp),intent(in),optional::tol_cov

      real(dp),allocatable::xm(:,:),xs(:,:),xn(:,:),xt(:,:),zbig(:,:,:)
      real(dp),allocatable::pm(:,:),ps(:,:),pn(:,:),pt(:,:),start(:)
      real(dp),allocatable::sigma(:,:,:),prec(:,:,:),snew(:,:),effects(:,:,:),postcov(:,:),pinv(:,:)
      integer,allocatable::levels(:)
      type(gamlss_control_t)::ctl
      integer::n,np,q,ng,j,it,nouter,istat,pfix(4),a,b,g,lo,hi
      real(dp)::crit,delta,maxdelta
      logical::act(4),corr(4)

      n=size(y);np=family_npar(family);q=size(z_random,2)
      if(n<=0.or.np<1.or.np>4.or.q<1.or.size(x_mu,1)/=n.or.size(group)/=n)then
         result%status=1;return
      end if
      if(size(z_random,1)/=n.or.size(z_random,3)<np)then;result%status=2;return;end if
      act=.false.;act(1)=.true.;if(present(active_parameters))act=active_parameters
      if(any(act(np+1:4)).or..not.any(act(1:np)))then;result%status=3;return;end if
      corr=.true.;if(present(correlated_within))corr=correlated_within
      call grouped_levels(group,levels,istat)
      if(istat/=0)then;result%status=4;return;end if
      ng=size(levels);result%levels=levels
      allocate(zbig(n,ng*q,4));zbig=0.0_dp
      do j=1,np
         if(act(j))call grouped_random_design_known(group,levels,z_random(:,:,j),zbig(:,:,j))
      end do
      call default_design(n,x_mu,x_sigma,x_nu,x_tau,xm,xs,xn,xt)
      pfix=[size(xm,2),size(xs,2),size(xn,2),size(xt,2)]
      call prepare_design(xm,zbig(:,:,1),act(1),pm)
      call prepare_design(xs,zbig(:,:,2),act(2),ps)
      call prepare_design(xn,zbig(:,:,3),act(3),pn)
      call prepare_design(xt,zbig(:,:,4),act(4),pt)

      allocate(sigma(q,q,4),prec(q,q,4));sigma=0.0_dp;prec=0.0_dp
      do j=1,np
         if(.not.act(j))cycle
         do a=1,q;sigma(a,a,j)=1.0_dp;prec(a,a,j)=1.0_dp;end do
      end do
      ctl=gamlss_control_t();if(present(control))ctl=control
      ctl%estimate_lambda_mu=.false.;ctl%estimate_lambda_sigma=.false.
      ctl%estimate_lambda_nu=.false.;ctl%estimate_lambda_tau=.false.
      nouter=20;if(present(max_outer))nouter=max(1,max_outer)
      crit=1.0e-5_dp;if(present(tol_cov))crit=max(0.0_dp,tol_cov)
      allocate(effects(ng,q,4));effects=0.0_dp

      do it=1,nouter
         call install_all_penalties(pfix,ng,q,act,prec,pm,ps,pn,pt)
         if(it==1)then
            call fit_gamlss_model(y,xm,family,result%model,x_sigma=xs,x_nu=xn,x_tau=xt,weights=weights, &
               offset_mu=offset_mu,offset_sigma=offset_sigma,offset_nu=offset_nu,offset_tau=offset_tau, &
               penalty_mu=pm,penalty_sigma=ps,penalty_nu=pn,penalty_tau=pt,lambda_mu=merge(1.0_dp,0.0_dp,act(1)), &
               lambda_sigma=merge(1.0_dp,0.0_dp,act(2)),lambda_nu=merge(1.0_dp,0.0_dp,act(3)), &
               lambda_tau=merge(1.0_dp,0.0_dp,act(4)),control=ctl)
         else
            call flatten_start(result%model,np,start)
            call fit_gamlss_model(y,xm,family,result%model,x_sigma=xs,x_nu=xn,x_tau=xt,weights=weights, &
               offset_mu=offset_mu,offset_sigma=offset_sigma,offset_nu=offset_nu,offset_tau=offset_tau, &
               penalty_mu=pm,penalty_sigma=ps,penalty_nu=pn,penalty_tau=pt,lambda_mu=merge(1.0_dp,0.0_dp,act(1)), &
               lambda_sigma=merge(1.0_dp,0.0_dp,act(2)),lambda_nu=merge(1.0_dp,0.0_dp,act(3)), &
               lambda_tau=merge(1.0_dp,0.0_dp,act(4)),start=start,control=ctl)
         end if
         if(result%model%status/=0)then;result%status=10+result%model%status;return;end if
         maxdelta=0.0_dp
         do j=1,np
            if(.not.act(j))cycle
            call extract_parameter_block(result%model,j,pfix(j),ng,q,effects(:,:,j),postcov,istat)
            if(istat/=0)then;result%status=5;return;end if
            allocate(snew(q,q));snew=0.0_dp
            do g=1,ng
               do a=1,q
                  do b=1,q
                     lo=(g-1)*q+a;hi=(g-1)*q+b
                     snew(a,b)=snew(a,b)+effects(g,a,j)*effects(g,b,j)+postcov(lo,hi)
                  end do
               end do
            end do
            snew=snew/real(ng,dp)
            if(.not.corr(j))call diagonalize(snew)
            call stabilize_covariance(snew,1.0e-8_dp)
            delta=maxval(abs(snew-sigma(:,:,j))/(1.0_dp+abs(sigma(:,:,j))))
            maxdelta=max(maxdelta,delta)
            sigma(:,:,j)=0.5_dp*sigma(:,:,j)+0.5_dp*snew
            if(.not.corr(j))call diagonalize(sigma(:,:,j))
            call invert_matrix(sigma(:,:,j),pinv,istat)
            if(istat==0)prec(:,:,j)=pinv
            if(allocated(pinv))deallocate(pinv)
            deallocate(snew,postcov)
            if(istat/=0)then;result%status=6;return;end if
         end do
         result%iterations=it
         if(maxdelta<crit)then;result%converged=.true.;exit;end if
      end do
      result%effects=effects;result%covariance=sigma;result%precision=prec
      result%active=act;result%correlated_within=corr;result%status=0
   end subroutine fit_gamlss_multi_random_effects

   subroutine grouped_levels(group,levels,status)
      integer,intent(in)::group(:)
      integer,allocatable,intent(out)::levels(:)
      integer,intent(out)::status
      integer,allocatable::tmp(:)
      integer::i,nlev
      if(size(group)==0)then;allocate(levels(0));status=1;return;end if
      allocate(tmp(size(group)));nlev=0
      do i=1,size(group)
         if(nlev==0)then;nlev=1;tmp(1)=group(i)
         else if(.not.any(tmp(1:nlev)==group(i)))then;nlev=nlev+1;tmp(nlev)=group(i);end if
      end do
      allocate(levels(nlev));levels=tmp(1:nlev);status=0
   end subroutine grouped_levels

   subroutine grouped_random_design_known(group,levels,z,zbig)
      integer,intent(in)::group(:),levels(:)
      real(dp),intent(in)::z(:,:)
      real(dp),intent(out)::zbig(:,:)
      integer::i,g,q
      q=size(z,2);zbig=0.0_dp
      do i=1,size(group)
         do g=1,size(levels)
            if(group(i)==levels(g))then;zbig(i,(g-1)*q+1:g*q)=z(i,:);exit;end if
         end do
      end do
   end subroutine grouped_random_design_known

   subroutine default_design(n,xmu,xs_in,xn_in,xt_in,xm,xs,xn,xt)
      integer,intent(in)::n
      real(dp),intent(in)::xmu(:,:)
      real(dp),intent(in),optional::xs_in(:,:),xn_in(:,:),xt_in(:,:)
      real(dp),allocatable,intent(out)::xm(:,:),xs(:,:),xn(:,:),xt(:,:)
      xm=xmu
      if(present(xs_in))then;xs=xs_in;else;allocate(xs(n,1));xs=1.0_dp;end if
      if(present(xn_in))then;xn=xn_in;else;allocate(xn(n,1));xn=1.0_dp;end if
      if(present(xt_in))then;xt=xt_in;else;allocate(xt(n,1));xt=1.0_dp;end if
   end subroutine default_design

   subroutine prepare_design(x,z,active,penalty)
      real(dp),allocatable,intent(inout)::x(:,:)
      real(dp),intent(in)::z(:,:)
      logical,intent(in)::active
      real(dp),allocatable,intent(out)::penalty(:,:)
      real(dp),allocatable::tmp(:,:)
      integer::p
      p=size(x,2)
      if(active)then
         allocate(tmp(size(x,1),p+size(z,2)));tmp(:,1:p)=x;tmp(:,p+1:)=z;call move_alloc(tmp,x)
      end if
      allocate(penalty(size(x,2),size(x,2)));penalty=0.0_dp
   end subroutine prepare_design

   subroutine install_all_penalties(pfix,ng,q,act,prec,pm,ps,pn,pt)
      integer,intent(in)::pfix(4),ng,q
      logical,intent(in)::act(4)
      real(dp),intent(in)::prec(:,:,:)
      real(dp),intent(inout)::pm(:,:),ps(:,:),pn(:,:),pt(:,:)
      integer::j,g,lo,hi
      pm=0.0_dp;ps=0.0_dp;pn=0.0_dp;pt=0.0_dp
      do j=1,4
         if(.not.act(j))cycle
         do g=1,ng
            lo=pfix(j)+(g-1)*q+1;hi=lo+q-1
            select case(j)
            case(1);pm(lo:hi,lo:hi)=prec(:,:,j)
            case(2);ps(lo:hi,lo:hi)=prec(:,:,j)
            case(3);pn(lo:hi,lo:hi)=prec(:,:,j)
            case(4);pt(lo:hi,lo:hi)=prec(:,:,j)
            end select
         end do
      end do
   end subroutine install_all_penalties

   subroutine extract_parameter_block(model,j,pfix,ng,q,effects,postcov,status)
      type(gamlss_result_t),intent(in)::model
      integer,intent(in)::j,pfix,ng,q
      real(dp),intent(out)::effects(:,:)
      real(dp),allocatable,intent(out)::postcov(:,:)
      integer,intent(out)::status
      real(dp),allocatable::coef(:),cov(:,:)
      integer::g,lo,hi
      status=0
      select case(j)
      case(1);coef=model%mu%coefficients;cov=model%mu%covariance
      case(2);coef=model%sigma%coefficients;cov=model%sigma%covariance
      case(3);coef=model%nu%coefficients;cov=model%nu%covariance
      case(4);coef=model%tau%coefficients;cov=model%tau%covariance
      case default;status=1;return
      end select
      if(size(coef)<pfix+ng*q.or.size(cov,1)<pfix+ng*q)then;status=2;return;end if
      allocate(postcov(ng*q,ng*q));postcov=0.0_dp
      do g=1,ng
         lo=pfix+(g-1)*q+1;hi=lo+q-1
         effects(g,:)=coef(lo:hi)
         postcov((g-1)*q+1:g*q,(g-1)*q+1:g*q)=cov(lo:hi,lo:hi)
      end do
   end subroutine extract_parameter_block

   subroutine flatten_start(model,np,start)
      type(gamlss_result_t),intent(in)::model
      integer,intent(in)::np
      real(dp),allocatable,intent(out)::start(:)
      integer::p(4),j,pos,total
      p=0;p(1)=size(model%mu%coefficients)
      if(np>=2)p(2)=size(model%sigma%coefficients)
      if(np>=3)p(3)=size(model%nu%coefficients)
      if(np>=4)p(4)=size(model%tau%coefficients)
      total=sum(p);allocate(start(total));pos=0
      do j=1,np
         select case(j)
         case(1);start(pos+1:pos+p(j))=model%mu%coefficients
         case(2);start(pos+1:pos+p(j))=model%sigma%coefficients
         case(3);start(pos+1:pos+p(j))=model%nu%coefficients
         case(4);start(pos+1:pos+p(j))=model%tau%coefficients
         end select
         pos=pos+p(j)
      end do
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
      a=0.5_dp*(a+transpose(a));do i=1,size(a,1);a(i,i)=max(a(i,i),eps);end do
   end subroutine stabilize_covariance
end module gamlss_multi_random_v05
